#!/usr/bin/env python3
#claridge_deduper_v2.py

import argparse
import re
import sys
import textwrap

parser = argparse.ArgumentParser(description="Removes PCR duplicates from a SAM file of sorted (using samtools) and uniquely mapped reads (single- and paired-end). For each set of duplicates, one duplicate is printed to an output file. Given a list of known UMIs, alignments with unexpected UMIs are ignored. If the reads are paired, assumes reads are mapped in proper pairs", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f','--file', help='absolute path to SAM file to be deduped', required=True, type=str)
parser.add_argument('-p','--paired', help='if flag is set, indicates that reads are paired-end (default: single-end)', required=False, action='store_true', default=False)
parser.add_argument('-u','--umi', help='absolute path to list of UMIs (default: randomers)', required=False, type=str)
parser.add_argument('-q','--quality_filter', help=textwrap.dedent('''if flag is set, duplicates will be quality filtered (default: first encountered duplicate)
note that for paired-end reads, quality filter is dictated by forward read only
options:
perbase   select duplicate with the highest per-base average quality
mapq      select duplicate with the highest MAPQ value'''), required=False, default=False)

args = parser.parse_args()



##########################################################################################
##### Define Higher Order Functions ######################################################
##########################################################################################

def strand_checker(flag):
    '''Takes the bitwise flag and checks it for strandedness.
    Assumes read is mapped, otherwise returnes None.
    Assumes data are single-end. Returns "+" or "-", depending on strand.'''
    
    # Check read is mapped
    if ((flag & 4) == 4):
        raise NameError("ERROR: Exiting program. Read is unmapped.")
    
    # Assume positive strand
    strand = "+"
    if ((flag & 16) == 16):
        # Changes strand if bit 16 is set
        strand = "-"
    return strand

def paired_read_checker(flag):
    '''Takes the bitwise flag.
    Checks if reads are paired-end and halts program if reads are not.
    Checks if reads are mapped in proper pairs and halts program if not.
    Checks if the read is forward or reverse and returns "forward" or "reverse".'''
    
    # Check if reads are paired
    if ((flag & 1) != 1):
        raise ValueError("ERROR: Exiting program. This is not paired-end data.")
    
    # Check if reads are mapped in proper pairs
    if ((flag & 2) != 2):
        raise ValueError("ERROR: Exiting program. Encountered reads that were not mapped in proper pairs.")
    
    # Check forward or reverse
    # Flags drived from Sequence Alignment/Map Format Specification (21 Aug 2017)
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    
    if ((flag & 64) == 64):
        read = "forward"
    if ((flag & 128) == 128):
        read = "reverse"
    return read

def POS_correct(cigar_string, POS):
    '''Takes the CIGAR string and corrects the start position in the SAM file (POS)
    according to any soft clipping at the beginning of the CIGAR string.
    Returns the corrected POS value.'''
    
    # Find soft clippling in begining of cigar_string
    search = re.search(r"^\d+S", cigar_string)
    
    # If soft clipping found
    if search:
        soft_clip = int(search.group(0)[:-1])
        corrected_POS = POS - soft_clip
    
    # If soft clipping not found
    else:
        corrected_POS = POS
        # Corrected starting position could be the same as the starting position
    return corrected_POS

def paired_corrections(cigar_string_f, cigar_string_r, POS_f, TLEN_f):
    '''For paired-end reads.
    Calculates corrected_POS for the first read via POS_correct function.
    Adds the TLEN to the corrected_POS, yielding the putative end of the template.
    Takes the CIGAR string of the second read and corrects the TLEN if soft clipping occurred.
    Returns the corrected TLEN value.'''
    corrected_POS = POS_correct(cigar_string_f, POS_f)
    end_TLEN = corrected_POS + TLEN_f
	
	# Find soft clippling in begining of cigar_string2
    search = re.search(r"^\d+S", cigar_string_r)
    
    # If soft clipping found
    if search:
        soft_clip = int(search.group(0)[:-1])
        corrected_TLEN = end_TLEN + soft_clip
    
    # If soft clipping not found
    else:
        corrected_TLEN = end_TLEN
	
    return [corrected_POS, corrected_TLEN]

def perbase_qscore(quality):
    '''Calculates the average per-base quality score (assumes Phred33 encoding) 
    based on the QUAL string in the SAM entry.'''
    score_total = 0
    for char in quality:
        score_total += ord(char) - 33 # Convert ASCII score to phred score
    score = score_total / len(quality)
    return score

def count_lines(infile):
    '''Opens the input file and returns the number of lines in the file.'''
    with open(infile) as file:
        for i, line in enumerate(file):
            pass
    return i + 1



##########################################################################################
##### Set Up #############################################################################
##########################################################################################

### Initialize empty dictionary for tuples
# tuple_dict will contain tuples as keys and SAM entries with highest average per-base quality as values
# Format of single-end tuples: (umi/randomer, strand +/-, chromosome, new_POS)
# Format of paired-end tuples: (umi/randomer pair, strand +/-, chromosome, new_POS, new_TLEN)
tuple_dict = {}

### Initialize empty set for UMIs if flag is set
# Make UMI set if UMI flag is set, otherwise leave empty
if args.umi is not None:
    umi_set = set()
    print("Creating UMI set:")
    with open(args.umi, "r+") as umi_list:
        for umi in umi_list:
            umi = umi.strip("\n")
            print("\t" + umi)
            umi_set.add(umi)
    

### Initialze first_line variable used in chromosome check
first_line = True

### Create output filename
infile = args.file
filename = infile.split("/")[-1].split(".")[0]   # Isolate file name
deduped = "".join([filename, "_deduped"])        # Add "deduped"
outfile = infile.replace(filename, deduped)      # Replace in original filepath

### Get file linecount to be used in last line check
total_lines = count_lines(infile)

print("Opening %s.sam for reading and %s.sam for writing." % (filename, deduped))

if args.quality_filter == False:
	print("Quality filter: Default.")
elif args.quality_filter == "perbase":
	print("Quality filter: Average per-base quality score.") 
elif args.quality_filter == "mapq": 
	print("Quality filter: mapq.") 
	
	
	
##########################################################################################
##### Read Through SAM File ##############################################################
##########################################################################################

##### Single-End #########################################################################

if args.paired == False:
	print("Read type: Single-end.")
	with open(infile, "r+") as file, open(outfile, "w+") as out:
		linecount = 0
		while True:
			line = file.readline()
		
			### Last line check
			if line == "":
				if args.quality_filter == False:
					pass
				if args.quality_filter != False:
					for value in tuple_dict.values():
						out.write(value)
				break
		
			linecount += 1
		
			### Print out header lines
			if line.startswith("@") == True:
				out.write(line)
		
			### Process alignments
			elif line.startswith("@") == False:
				line_list = line.strip("\n").split("\t")
		
				### Progress report
				if linecount % 200000 == 0:
					print("Passed line " + str(linecount) + ".")
		  
				### Check bitwise flag
				flag = int(line_list[1])
				# Check strand
				strand = strand_checker(flag)

				### Check UMI
				umi = line_list[0].split(":")[-1]
				if args.umi is None:                  # If UMI flag was not set, assume randomers and do nothing
					pass
				elif args.umi is not None:            # If UMI flag was set
					if umi in umi_set:				  # If UMI is in the umi_set, do nothing
						pass
					else:
						if linecount == total_lines:      # If last line in file, print tuple_dict contents and break
							for value in tuple_dict.values():
								out.write(value)
							break
						continue                          # Do not write unexpexted alignment with unexpected UMI to file
		   
				### Check chromosome
				# If first line, set the chromosome variable
				if first_line == True:
					chromosome = line_list[2]
					first_line == False
				# If not the first line, run chromosome check for emptying tuple_dict
				elif first_line == False:
					prev_chromosome = chromosome
					chromosome = line_list[2]
					if chromosome != prev_chromosome:  # Inequality means that you have switched chromosomes
						for value in tuple_dict.values():
							out.write(value)           # Print tuple_dict contents to output file
						tuple_dict = {}                # Empty tuple_dict
		
				### Check POS and CIGAR string (correcting for soft clipping if necessary)
				POS = int(line_list[3])
				cigar_string = line_list[5]
				new_POS = POS_correct(cigar_string, POS)
			
				### Make tuple
				query_tuple = (umi, strand, chromosome, new_POS)				
			
				### Check tuple_dict, method depending on quality_filter flag
				if args.quality_filter == False:
					if query_tuple not in tuple_dict:
						tuple_dict[query_tuple] = 1
						out.write(line)
					elif query_tuple in tuple_dict:
						continue
				elif args.quality_filter == "perbase": 	
					# Calculate average per-base quality score across the alignment
					quality = line_list[10]
					score = perbase_qscore(quality)
					if query_tuple in tuple_dict:
						# Check existing score against new alignment's score, replace if new alignments's is higher
						existing_quality = tuple_dict[query_tuple].split("\t")[10]
						existing_score = perbase_qscore(existing_quality)
						if existing_score >= score:
							pass
						elif existing_score < score:
							tuple_dict[query_tuple] = line
					if query_tuple not in tuple_dict:
						# Add alignment if tuple not found in tuple_dict
						tuple_dict[query_tuple] = line
				elif args.quality_filter == "mapq": 	
					MAPQ = line_list[4]
					if query_tuple in tuple_dict:
						# Check existing MAPQ against new alignment's MAPQ, replace if new alignments's is higher
						existing_MAPQ = tuple_dict[query_tuple].split("\t")[4]
						if existing_MAPQ >= MAPQ:
							pass
						elif existing_MAPQ < MAPQ:
							tuple_dict[query_tuple] = line
					if query_tuple not in tuple_dict:
						# Add alignment if tuple not found in tuple_dict
						tuple_dict[query_tuple] = line
					
##### Paired-End #########################################################################

if args.paired == True:
	print("Read type: Paired-end.")
	with open(infile, "r+") as file, open(outfile, "w+") as out:
		linecount = 0
		header_count = 0
		singleton_dict = {}
		while True:
			line1 = file.readline()
			linecount += 1
			
			### Progress report
			if linecount % 200000 == 0:
				print("Passed line " + str(linecount) + ".")
				
			### Print out header lines
			if line1.startswith("@") == True:
				out.write(line1)
		
			### Process alignments
			elif line1.startswith("@") == False:
				### Last line check
				if line1 == "":
					if args.quality_filter == False:
						pass
					if args.quality_filter != False:
						for value in tuple_dict.values():
							out.write(value[0])
							out.write(value[1])
					break
				
				### Check if paired reads are adjacent in sorted SAM file
				line1_list = line1.strip("\n").split("\t")
				QNAME = line1_list[0]			
				if QNAME not in singleton_dict:
					singleton_dict[QNAME] = line1
					continue
				elif QNAME in singleton_dict:
					line2 = singleton_dict.pop(QNAME)	      # Remove singleton from dictionary and save value to new variable
					line2_list = line2.strip("\n").split("\t")
					
					### Check UMI
					umi = line1_list[0].split(":")[-1]
					if args.umi is None:                      # If UMI flag was not set, assume randomers and do nothing
						pass
					elif args.umi is not None:                # If UMI flag was set
						umi_a = umi.split("^")[0]
						umi_b = umi.split("^")[1]
						if umi_a in umi_set and umi_b in umi_set:				  # If both UMIs are in the umi_set, do nothing
							pass
						else:
							if linecount >= total_lines:      # If last line in file or an empty line, print tuple_dict contents and break
								for value in tuple_dict.values():
									out.write(value)
								break
							continue                          # Do not write unexpexted alignment with unexpected UMI to file
					
					### Check bitwise flags
					flag1 = int(line1_list[1])
					flag2 = int(line2_list[1])

					# Check read for line1 and line2
					read1 = paired_read_checker(flag1)
					read2 = paired_read_checker(flag2)
					# Reassign line and line_list variable names to f (forward) and r (reverse) based on read
					if read1 == "reverse" and read2 == "forward":
						line_f = line2
						line_r = line1
						flag_f = flag2
						line_f_list = line2_list
						line_r_list = line1_list
					elif read2 == "reverse" and read1 == "forward":
						line_f = line1
						line_r = line2
						flag_f = flag1
						line_f_list = line1_list
						line_r_list = line2_list
					# Check strand
					strand = strand_checker(flag_f)

					### Check chromosome
					# If first line, set the chromosome variable
					if first_line == True:
						chromosome = line_f_list[2]
						first_line == False
					# If not the first line, run chromosome check for emptying tuple_dict
					elif first_line == False:
						prev_chromosome = chromosome
						chromosome = line_f_list[2]
						if chromosome != prev_chromosome:  # Inequality means that you have switched chromosomes
							for value in tuple_dict.values():
								out.write(value)           # Print tuple_dict contents to output file
							tuple_dict = {}                # Empty tuple_dict
		
					### Check POS, TLEN, and CIGAR string (correcting for soft clipping if necessary)
					cigar_string_f = line_f_list[5]
					cigar_string_r = line_r_list[5]
					POS_f = int(line_f_list[3])
					TLEN_f = int(line_f_list[8])
					# Produce corrected_POS and corrected_TLEN with paired_corrections function
					corrections = paired_corrections(cigar_string_f, cigar_string_r, POS_f, TLEN_f)
					new_POS = corrections[0]
					new_TLEN = corrections[1]
					
					### Make tuple
					query_tuple = (umi, strand, chromosome, new_POS, new_TLEN)				
			
					### Check tuple_dict, method depending on quality_filter flag
					if args.quality_filter == False:
						if query_tuple not in tuple_dict:
							tuple_dict[query_tuple] = 1
							out.write(line_f)
							out.write(line_r)
						elif query_tuple in tuple_dict:
							continue
					elif args.quality_filter == "perbase":  	
						# Calculate average per-base quality score across the forward read
						quality = line_f_list[10]
						score = perbase_qscore(quality)
						if query_tuple in tuple_dict:
							# Check existing score against new alignment's score, replace if new alignments's is higher
							existing_quality = tuple_dict[query_tuple][0].split("\t")[10]
							existing_score = perbase_qscore(existing_quality)
							if existing_score >= score:
								pass
							elif existing_score < score:
								tuple_dict[query_tuple] = [line_f, line_r]
						if query_tuple not in tuple_dict:
							# Add reads to tuple_dict if tuple not found in tuple_dict
							tuple_dict[query_tuple] = [line_f, line_r]
					elif args.quality_filter == "mapq":	
						MAPQ = line_f_list[4]
						if query_tuple in tuple_dict:
							# Check existing MAPQ against new alignment's MAPQ, replace if new alignments's is higher
							existing_MAPQ = tuple_dict[query_tuple][0].split("\t")[4]
							if existing_MAPQ >= MAPQ:
								pass
							elif existing_MAPQ < MAPQ:
								tuple_dict[query_tuple] = [line_f, line_r]
						if query_tuple not in tuple_dict:
							# Add alignment if tuple not found in tuple_dict
							tuple_dict[query_tuple] = [line_f, line_r]
							
print("Finished.")
sys.exit()

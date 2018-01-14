#!/usr/bin/env python3

# claridge_deduper.py
# Sally Claridge
# Due 15 January 2018

import argparse
import re
import sys
import textwrap

parser = argparse.ArgumentParser(description="Removes PCR duplicates from a SAM file (with .sam extension) of samtools sorted and uniquely mapped reads (single- and paired-end). For each set of duplicates, one duplicate is printed to an output file. Given a list of known UMIs, alignments with unexpected UMIs are ignored. If the reads are paired, assumes reads are mapped in proper pairs", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f','--file', help='absolute path to SAM file to be deduped', required=True, type=str)
parser.add_argument('-p','--paired', help='if flag is set, indicates that reads are paired-end (default: single-end)', required=False, action='store_true', default=False)
parser.add_argument('-u','--umi', help='absolute path to list of UMIs (default: randomers)', required=False, type=str)
parser.add_argument('-n','--name', help='name output <FILENAME>_deduped.sam (default: <FILENAME>.sam_deduped)', default=False, action='store_true', required=False)
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
    
    # Change strand if bit 16 is set
    if ((flag & 16) == 16):
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
    
    # Check for forward or reverse
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
    
    # If soft clipping found, correct POS to reflect soft-clipped bases
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
    
    # Correct POS of forward read
    corrected_POS = POS_correct(cigar_string_f, POS_f)
    
    # Add TLEN to corrected POS
    end_TLEN = corrected_POS + TLEN_f
	
	# Find soft clippling in begining of cigar_string2
    search = re.search(r"^\d+S", cigar_string_r)
    
    # If soft clipping found, correct TLEN to reflect soft-clipped bases
    if search:
        soft_clip = int(search.group(0)[:-1])
        corrected_TLEN = end_TLEN + soft_clip
    
    # If soft clipping not found
    else:
        corrected_TLEN = end_TLEN
	
    return [corrected_POS, corrected_TLEN]

def perbase_qscore(quality):
    '''Takes in the QUAL string in the SAM alignment.
    Calculates the string's average per-base quality score (assumes Phred33 encoding).'''
    
    # Initialize score counter
    score_total = 0
    
    # Convert each ASCII score to phred score and add to score counter
    for char in quality:
        score_total += ord(char) - 33
    
    # Calculate average score per base
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

### Initialize empty set of UMIs if umi flag is set
# Print out UMIs to stdout for reference
if args.umi is not None:
    umi_set = set()
    print("Creating UMI set:")
    with open(args.umi, "r+") as umi_list:
        for umi in umi_list:
            umi = umi.strip("\n").upper()		   	 # Make uppercase in case lowercase UMIs were provided
            print("\t" + umi)
            umi_set.add(umi)
    

### Initialze first_line variable used in chromosome check
first_line = True

### Create output filename
infile = args.file

# If name flag is set (format "<FILENAME>_deduped.sam")
if args.name == True:
	filename = infile.split("/")[-1].split(".")[0]   # Isolate file name without .sam
	deduped = "".join([filename, "_deduped"])        # Add "deduped" to filename
	outfile = infile.replace(filename, deduped)      # Replace in original filepath
	print("Opening %s.sam for reading and %s.sam for writing." % (filename, deduped))

# If name flag is not set (format "<FILENAME>.sam_deduped")
elif args.name == False:
	filename = infile.split("/")[-1]				 # Isolate faile name with .sam
	deduped = "".join([filename, "_deduped"])		 # Add "deduped" to filename
	outfile = "".join([infile, "_deduped"])			 # Add "deduped" to original filepath
	print("Opening %s for reading and %s for writing." % (filename, deduped))

### Get file linecount to be used in last line check
total_lines = count_lines(infile)

### Print quality filter specification
if args.quality_filter == False:
	print("Quality filter: Default, first encountered duplicate.")
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
		head_count = 0
		output_count = 0
		duplicate_count = 0
		umi_n_count = 0
		unexpected_umi_count = 0
		
		while True:
			line = file.readline()
		
			### Last line check
			if line == "":
				if args.quality_filter == False:			# Without quality_filter, do nothing
					print("Wrote the chromosome %s set of alignments to output file." % chromosome)
					pass
				if args.quality_filter != False:			# With quality_filter, print tuple_dict contents before breaking
					for value in tuple_dict.values():
						out.write(value)
						output_count += 1
					print("Wrote the chromosome %s set of alignments to output file." % chromosome)
				break
			
			# Increment line counter for each line encountered
			linecount += 1
		
			### Print out header lines
			if line.startswith("@") == True:
				out.write(line)
				head_count += 1
		
			### Process alignments
			elif line.startswith("@") == False:
				line_list = line.strip("\n").split("\t")	# Split each alignment by tabs
		
				### Progress report
				if linecount % 200000 == 0:
					print("Passed line " + str(linecount) + ".")
		  
				### Check bitwise flag (convert to int type)
				flag = int(line_list[1])
				# Check strand
				strand = strand_checker(flag)

				### Check UMI
				# Assumes UMI is at the end of QNAME
				umi = line_list[0].split(":")[-1]
				if args.umi is None:                  		# If UMI flag was not set, assume randomers
					if bool(re.search("N", umi)) == True:	# If randomer contains N, increment N counter, do nothing with this alignment and move to next line
						umi_n_count += 1
						continue
					elif bool(re.search("N", umi)) == False:
						pass								# Else pass and do nothing
				elif args.umi is not None:            		# If UMI flag was set, check it against the umi_set
					if umi in umi_set:				  		# If UMI is in the umi_set, do nothing
						pass
					else:									# If UMI is not in the umi_set
						if bool(re.search("N", umi)) == True:	# If UMI contains N, increment N counter
							umi_n_count += 1
						if bool(re.search("N", umi)) == False:	# If umi does not contain N, increment unexpected counter
							unexpected_umi_count += 1
						if linecount >= total_lines:      	# If last line in file, print tuple_dict contents and break out of while loop
							for value in tuple_dict.values():
								out.write(value)
								output_count += 1
							print("Wrote the chromosome %s set of alignments to output file." % chromosome)
							break
						continue                          	# Do not write unexpected alignment with unexpected UMI to file
		   
				### Check chromosome
				if first_line:								# If first line, set the chromosome variable
					chromosome = line_list[2]
					first_line == False
				elif not first_line:						# If not the first line, run chromosome check for emptying tuple_dict
					prev_chromosome = copy.copy(chromosome)	# Reassign chromosome to new variable
					chromosome = line_list[2]				# Set chromosome variable for new line
					if chromosome != prev_chromosome:  		# Inequality means that you have switched chromosomes
						for value in tuple_dict.values():
							out.write(value)           		# Print tuple_dict contents to output file
							output_count += 1
						print("Wrote the chromosome %s set of alignments to output file." % prev_chromosome)
						tuple_dict = {}                		# Empty tuple_dict
					else:
						pass
		
				### Check POS and CIGAR string (correcting for soft clipping)
				POS = int(line_list[3])
				cigar_string = line_list[5]
				new_POS = POS_correct(cigar_string, POS)	# Correct for soft clipping
			
				### Make tuple identifier
				query_tuple = (umi, strand, chromosome, new_POS)				
			
				### Check tuple_dict, method depending on quality_filter flag
				# No quality filtering
				if args.quality_filter == False:			# Output first encountered duplicate
					if query_tuple not in tuple_dict:		# If tuple identifier has not been seen before
						tuple_dict[query_tuple] = 1			# Add this new tuple identifier to the tuple_dict
						out.write(line)						# Write
						output_count += 1
					elif query_tuple in tuple_dict:			# Go to next line if this tuple identifier has been seen before (duplicate)
						duplicate_count += 1
						continue
				
				# Per-base quality filtering
				elif args.quality_filter == "perbase": 		# Calculate average per-base quality score across the alignment
					quality = line_list[10]
					score = perbase_qscore(quality)
					if query_tuple in tuple_dict:
						duplicate_count += 1
						# Check existing score against new alignment's score, replace if new alignments's is higher
						existing_quality = tuple_dict[query_tuple].split("\t")[10]
						existing_score = perbase_qscore(existing_quality)
						if existing_score >= score:
							pass							# Do not replace
						elif existing_score < score:
							tuple_dict[query_tuple] = line	# Replace
					if query_tuple not in tuple_dict:		# Add alignment if tuple not found in tuple_dict
						tuple_dict[query_tuple] = line
				
				# MAPQ quality filtering
				elif args.quality_filter == "mapq": 	
					MAPQ = line_list[4]
					if query_tuple in tuple_dict:
						duplicate_count += 1
						# Check existing MAPQ against new alignment's MAPQ, replace if new alignments's is higher
						existing_MAPQ = tuple_dict[query_tuple].split("\t")[4]
						if existing_MAPQ >= MAPQ:
							pass							# Do nothing
						elif existing_MAPQ < MAPQ:
							tuple_dict[query_tuple] = line	# Replace
					if query_tuple not in tuple_dict:		# Add alignment if tuple not found in tuple_dict
						tuple_dict[query_tuple] = line
					
##### Paired-End #########################################################################

if args.paired == True:
	print("Read type: Paired-end.")
	with open(infile, "r+") as file, open(outfile, "w+") as out:
		linecount = 0
		head_count = 0
		output_count = 0
		duplicate_count = 0
		umi_n_count = 0
		unexpected_umi_count = 0
				
		### Initialize dictionary to contain unpaired reads based on QNAME
		# Used to check if paired reads are adjacent in sorted SAM file
		# Key is QNAME, value is whole SAM alignment
		singleton_dict = {}
		
		while True:
			line1 = file.readline()
			
			# Increment line counter for each line encountered
			linecount += 1
			
			### Progress report
			if linecount % 200000 == 0:
				print("Passed line " + str(linecount) + ".")

			### Print out header lines
			if line1.startswith("@") == True:
				out.write(line1)
				head_count += 1
		
			### Process alignments
			elif line1.startswith("@") == False:
				
				### Last line check
				if line1 == "":
					if args.quality_filter == False:			# Without quality_filter, do nothing
						print("Wrote the chromosome %s set of alignments to output file." % chromosome)
						pass
					if args.quality_filter != False:			# With quality_filter, print tuple_dict contents before breaking
						for value in tuple_dict.values():
							out.write(value[0])					# Forward read
							out.write(value[1])					# Reverse read
							output_count += 2
						print("Wrote the chromosome %s set of alignments to output file." % chromosome)
					break
				
				### Check if data is paired by running paired_read_checker on flag
				# Breaks if data is not paired-end
				line1_list = line1.strip("\n").split("\t")		# Split alignment by tabs
				flag1 = int(line1_list[1])
				read = paired_read_checker(flag1)
				
				### Check if paired reads are adjacent in sorted SAM file
				QNAME = line1_list[0]							# Save QNAME, which should pair with only one other read in the SAM file	
				if QNAME not in singleton_dict:					# Save read to singleton_dict if it's pair hasn't been seen yet
					singleton_dict[QNAME] = line1
					continue
				elif QNAME in singleton_dict:
					line2 = singleton_dict.pop(QNAME)	      	# Remove singleton from dictionary and save value to new variable
					line2_list = line2.strip("\n").split("\t")	# Split alignment by tabs
					
					### Check UMI
					umi = line1_list[0].split(":")[-1]
					if args.umi is None:                  		# If UMI flag was not set, assume randomers
						if bool(re.search("N", umi)) == True:	# If either randomer in pair contains N, do nothing with this alignment and move to next line
							umi_n_count += 1
							continue
						elif bool(re.search("N", umi)) == False:
							pass								# Else pass and do nothing
					elif args.umi is not None:                	# If UMI flag was set, isolate both UMIs
						umi_a = umi.split("^")[0]
						umi_b = umi.split("^")[1]
						if umi_a in umi_set and umi_b in umi_set:
							pass								# If both UMIs are in the umi_set, do nothing
						else:
							if bool(re.search("N", umi)) == True:	# If either UMI in pair contains N, incremenet N counter
								umi_n_count += 1
							if bool(re.search("N", umi)) == False:	# If either UMI in pair does not contain N, incremenet unexpected counter
								unexpected_umi_count += 1
							if linecount >= total_lines:      	# If last line in file or an empty line, print tuple_dict contents and break
								for value in tuple_dict.values():
									out.write(value[0])			# Forward read
									out.write(value[1])			# Reverse read
									output_count += 2
								print("Wrote the chromosome %s set of alignments to output file." % chromosome)
								break
							continue                          	# Do not write unexpexted alignment with unexpected UMI to file
					
					### Check bitwise flags
					flag1 = int(line1_list[1])
					flag2 = int(line2_list[1])

					# Check read for line1 and line2 (forward or reverse) using paired_read_checker
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
					if first_line:								# If first line, set the chromosome variable
						chromosome = line_f_list[2]
						first_line == False
					elif not first_line:						# If not the first line, run chromosome check for emptying tuple_dict
						prev_chromosome = copy.copy(chromosome)	# Reassign chromosome to new variable
						chromosome = line_f_list[2]				# Set chromosome variable for new line
						if chromosome != prev_chromosome:  		# Inequality means that you have switched chromosomes
							for value in tuple_dict.values():	# Print tuple_dict contents to output file
								out.write(value[0])           	# Forward read
								out.write(value[1])				# Reverse read
								output_count += 2
							print("Wrote the chromosome %s set of alignments to output file." % prev_chromosome)
							tuple_dict = {}                		# Empty tuple_dict
						else:
							pass
		
					### Check POS, TLEN, and CIGAR string (correcting for soft clipping if necessary)
					cigar_string_f = line_f_list[5]
					cigar_string_r = line_r_list[5]
					POS_f = int(line_f_list[3])
					TLEN_f = int(line_f_list[8])
					
					# Produce corrected_POS and corrected_TLEN with paired_corrections function
					corrections = paired_corrections(cigar_string_f, cigar_string_r, POS_f, TLEN_f)
					new_POS = corrections[0]
					new_TLEN = corrections[1]
					
					### Make tuple identifier
					query_tuple = (umi, strand, chromosome, new_POS, new_TLEN)				
			
					### Check tuple_dict, method depending on quality_filter flag
					# No quality filtering
					if args.quality_filter == False:			# Output first encountered duplicate pair
						if query_tuple not in tuple_dict:		# If tuple identifier has not been seen before
							tuple_dict[query_tuple] = 1			# Add this new tuple identifier to the tuple_dict
							out.write(line_f)					# Write forward read
							out.write(line_r)					# Write reverse read
							output_count += 2
						elif query_tuple in tuple_dict:			# Go to next line if this tuple identifier has been seen before (duplicate)
							duplicate_count += 1
							continue
					
					# Per-base quality filtering
					elif args.quality_filter == "perbase":  	# Calculate average per-base quality score across the forward read
						quality = line_f_list[10]
						score = perbase_qscore(quality)
						if query_tuple in tuple_dict:
							duplicate_count += 1
							# Check existing score against new alignment's score, replace if new alignments's is higher
							existing_quality = tuple_dict[query_tuple][0].split("\t")[10]
							existing_score = perbase_qscore(existing_quality)
							if existing_score >= score:			# Do nothing
								pass
							elif existing_score < score:		# Replace
								tuple_dict[query_tuple] = [line_f, line_r]
						if query_tuple not in tuple_dict:
							# Add reads to tuple_dict if tuple not found in tuple_dict
							tuple_dict[query_tuple] = [line_f, line_r]
					
					# MAPQ quality filtering
					elif args.quality_filter == "mapq":	
						MAPQ = line_f_list[4]
						if query_tuple in tuple_dict:
							duplicate_count += 1
							# Check existing MAPQ against new alignment's MAPQ, replace if new alignments's is higher
							existing_MAPQ = tuple_dict[query_tuple][0].split("\t")[4]
							if existing_MAPQ >= MAPQ:			# Do nothing
								pass
							elif existing_MAPQ < MAPQ:			# Replace
								tuple_dict[query_tuple] = [line_f, line_r]
						if query_tuple not in tuple_dict:		# Add alignment if tuple not found in tuple_dict
							tuple_dict[query_tuple] = [line_f, line_r]



##########################################################################################
##### Closing Remarks ####################################################################
##########################################################################################		
					
print("Output file counts:")
print("\t" + "Total lines: " + str(head_count + output_count))
print("\t" + "Header lines: " + str(head_count))
if args.paired == False:
	print("\t" + "Alignments: " + str(output_count))
	print("\t" + "Duplicates removed: " + str(duplicate_count))
	if args.umi is None:
		print("\t" + "Randomers containing 1+ N: " + str(umi_n_count))
	if args.umi is not None:
		print("\t" + "UMIs containing 1+ N: " + str(umi_n_count))
		print("\t" + "Unexpected UMIs: " + str(unexpected_umi_count))
elif args.paired == True:
	print("\t" + "Alignments (pairs): " + str(output_count / 2))
	print("\t" + "Duplicates removed (pairs): " + str(duplicate_count))
	if args.umi is None:
		print("\t" + "Pairs of randomers containing 1+ N: " + str(umi_n_count))
	if args.umi is not None:
		print("\t" + "Pairs of UMIs containing 1+ N: " + str(umi_n_count))
		print("\t" + "Pairs containing unexpected UMIs: " + str(unexpected_umi_count))
	

print("Finished.")
sys.exit()

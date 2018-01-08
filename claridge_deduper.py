#!/usr/bin/env python
#claridge_deduper.py

import argparse
import re
import sys

parser = argparse.ArgumentParser(description="Removes PCR duplicates from a SAM file of uniquely mapped, single-end reads. For each set of duplicates, the duplicate with the highest per-base average quality is printed to an output file. Given a list of known UMIs, alignments with unexpected UMIs are ignored.")
parser.add_argument('-f','--file', help='absolute path to SAM file to be deduped.', required=True, type=str)
parser.add_argument('-p','--paired', help='reads are paired-end.', required=False, action='store_true', default=False)
parser.add_argument('-u','--umi', help='absolute path to list of UMIs (default: randomers were used).', required=False, type=str)
args = parser.parse_args()

# if args.paired:
#     print("Exiting program. No paired-end functionality at this time.")
#     sys.exit()


#####################################################################################
##### Define Higher Order Functions #################################################
#####################################################################################

def strand_checker(flag):
    '''Takes the bitwise flag and checks it for strandedness. Assumes read is mapped, otherwise returnes None.
    Assumes data are single-end. Returns "+" or "-", depending on strand.'''
    # Check read is mapped
    if ((flag & 4) == 4):
        raise NameError("Exiting program. Read is unmapped.")
        return None
    # Assume positive strand
    strand = "+"
    if ((flag & 16) == 16):
        # Changes strand if bit 16 is set
        strand = "-"
    return strand

def read_checker(flag):
    '''Takes the bitwise flag and checks whether the reads are paired-end or not.
    Halts program if reads are unpaired. Returns "forward" or "reverse", depending on the read.'''
    # Check if reads are paired
    if ((flag & 1) != 1):
        raise NameError("Exiting program. This is not paried-end data.")
        return None
    if ((flag & 40) == 40):
        read = "forward"
    if ((flag & 80) == 80):
        read = "reverse"
    return read

def POS_correct(cigar_string, POS):
    '''Takes the CIGAR string and corrects the start position in the SAM file (POS)
    according to any soft clipping at the beginning of the CIGAR string. Returns the corrected POS value.'''
    # Find soft clippling in begining of cigar string
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



#####################################################################################
##### Set Up ########################################################################
#####################################################################################

### Initialize empty dictionaries for tuples and dictionaries
tuple_dict = {}
umi_dict = {}

### Make UMI dictionary if UMI flag is set
if args.umi is not None:
    with open(args.umi) as umi_list:
        for umi in umi_list:
            umi = umi.strip("\n")
            umi_dict[umi] = 0
    print("Created UMI dictionary.", flush = True)

### Initialze chromosome and read variables
# For use in check to empty tuple_dict
chromosome = 1
# read will be a non-factor in single-end data, but will be edited if paired flag set
read = "na"

### Create output filename
infile = "/Users/sally_claridge/Desktop/deduper-sclaridg/test_input.sam"
filename = infile.split("/")[-1].split(".")[0]   # Isolate file name
deduped = "".join([filename, "_deduped"])        # Add "deduped"
outfile = infile.replace(filename, deduped)      # Replace in original filepath

### Get file linecount
total_lines = count_lines(infile)



#####################################################################################
##### Read Through SAM File #########################################################
#####################################################################################

print("Opening SAM files for reading and writing.", flush = True)

with open(infile, "r+") as file, open(outfile, "w+") as out:
    linecount = 0
    for line in file:
        linecount += 1
        if line.startswith("@") == True:
            pass
        elif line.startswith("@") == False:
            line_list = line.strip("\n").split("\t")
            
            ### Progress report
            if linecount % 200000 == 0:
                print("Passed line" + str(linecount) + ".", flush = True)
            
            ### Check UMI
            umi = line_list[0].split(":")[-1]
            if args.umi is None:                # If UMI flag was not set, assume randomers
                if umi in umi_dict:
                    umi_dict[umi] +=1
                else:
                    umi_dict[umi] = 1           # Still use umi_dict, but store randomers instead of UMIs
            elif args.umi is not None:          # If UMI flag was set
                if umi in umi_dict:
                    umi_dict[umi] +=1
            else:
                if linecount == total_lines:
                    for value in tuple_dict.values():
                        out.write(value)
                continue
                
            ### Check bitwise flag
            flag = int(line_list[1])
            # Check strand
            strand = strand_checker(flag)
            # Check read (if paired flag set)
            if args.paired:
                read = read_checker(flag)
                
            ### Check chromosome
            prev_chromosome = chromosome
            chromosome = line_list[2]
            if chromosome != prev_chromosome:
                for value in tuple_dict.values():
                    out.write(value)
                tuple_dict = {}
            
            ### Check POS and CIGAR string
            POS = int(line_list[3])
            cigar_string = line_list[5]
            new_POS = POS_correct(cigar_string, POS)
            
            ### Calculate average quality score across the alignment
            quality = line_list[10]
            score = perbase_qscore(quality)
            
            ### Check tuple_dict
            query_tuple = (umi, strand, chromosome, new_POS, read)
            if query_tuple in tuple_dict:
                existing_quality = tuple_dict[query_tuple].split("\t")[10]
                existing_score = perbase_qscore(existing_quality)
                if existing_score >= score:
                    pass
                elif existing_score < score:
                    tuple_dict[query_tuple] = line
            if query_tuple not in tuple_dict:
                tuple_dict[query_tuple] = line
                
            ### Last line check
            if linecount == total_lines:
                for value in tuple_dict.values():
                    out.write(value)

### Print closing remarks                    
print("\n")
print("Molecular_Identifier\tCount", flush = True)
for k, v in umi_dict.items():
    print(str(k) + "\t" + str(v), flush = True)

print("\n") 
print("Finished.", flush = True)
sys.exit()

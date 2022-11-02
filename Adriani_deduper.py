#!/usr/bin/env python

import re
import argparse
import bioinfo

def get_args():
    parser = argparse.ArgumentParser(description="This code will take a sorted, uniquely mapped SAM file with a known set of UMIs and return a deduplicated SAM file when the known UMI, Chromosome, Strand, and corrected Start position are all the same.")
    parser.add_argument("-f",  "--file", help="designates absolute file path to sorted sam file", required = True)
    parser.add_argument("-o", "--outfile", help = "designates absolute file path to sorted sam file", required = True)
    parser.add_argument("-u", "--umi", help = "designates file containing the list of UMIs", required = True)
    #parser.add_argument("-h", "--help", help = "This code takes uniquely mapped, samtools-sorted SAM file and deduplicates by known UMIs")
    return parser.parse_args()
	
args = get_args()

UMIset = set() 
#Create a set of the known UMIs
with open(args.umi, "r") as UMIS:
    for umi in UMIS: 
        umi = umi.strip()
        UMIset.add(umi)

requirements = dict()
header = 0
extras = 0
false_umi = 0
unique = 0
chromosome = dict()

SAM = open(args.file,"r")
OUT = open(args.outfile,"w")

for line in SAM: 
    if line.startswith("@")==True:
        OUT.write(line)
        header += 1
    else:  # a non starter line, does not start with @. 
        columns = line.split() #split by the tabs
        umi = bioinfo.get_UMI(columns[0]) #get the UMI out of the header
        if umi in UMIset: #if its in our known list of UMIS
            POS = columns[3] #save the position 
            cigar = columns[5] #save the cigar string
            strand = bioinfo.bitwise_strand_check(columns[1]) #get a "+" or "-" value for the strand, to put into correct_pos
            CHROM = columns[2] #get the chromsome
            POS = bioinfo.correct_pos(cigar,strand,POS) #correct it! depending on strand and cigar string, adjust position
            if (umi,POS,strand,CHROM) in requirements: #if this is already seen, increase count and don't write out
                requirements[(umi,POS,strand,CHROM)]+= 1
                extras +=1
            if (umi,POS,strand,CHROM) not in requirements: #if it's not in it, add it and write to the file. 
                requirements[(umi,POS,strand,CHROM)]= 1
                unique += 1
                if CHROM in chromosome:
                    chromosome[CHROM] += 1
                else :
                    chromosome[CHROM] = 1
                OUT.write(line)
        else: 
            false_umi += 1


with open("requirements.txt", "w") as r:
    r.write("Total Unique Reads found: "+str(unique)+"\n")
    r.write("Total Duplicates found: "+str(extras)+"\n")
    r.write("Header lines found: "+str(header)+"\n")
    r.write("Wrong UMIs found: "+str(false_umi)+"\n")
    r.write("Chromosomes seen in file"+"\n")
    for key, value in chromosome.items():
        count = f"{key}\t{value}\n"
        r.write(count)








# #def correct_pos(cigar,strand,position):
#     '''Takes the cigar string (as a string) and the position from the SAM file, and strand, and corrects for starting position'''
#     cig = str(cigar)
#     pos = int(position)
#     if strand == "+":
#         if "S" in cig:
#             cig = cig.split("S") 
#             part = bool(re.search(r"^[\d']", cig[1]))
#             if part==True:
#                 soft = int(cig[0])
#                 start_position = position + soft
#             else: 
#                 start_position = position
#         else: 
#             start_position = position
#     if strand == "-":
#         ## if there are letters?
#         if "M" in cig: 
#             cuts = cig.split("M")
#             position = 
#         if "D" in cig: 
#             position = 
#         if cig.endswith("S")==True: #this section works 
#             soft = cig.split("S")
#             soft = re.split('[\D]',soft[-2]) 
#             soft = int(soft[-1])
#             start_position = position + soft
#         # else: 
#         #     if "[A-Z]" in cig: 
#         #         cutting = re.split("\d",cig)

#         #     soft = re.search(r'M', cig) #find M in the cigar string. Split by soft?

        #if cig.endswith("S")==True:
            #soft = cig.split("S")
            # soft = re.split(r'\D+',soft) #then it should cut where the number ends
            # soft = cig.split[-2]
            # if "[A-Z]" in soft: 
            #     soft = cig.split("[A-Z]")

            # soft = re.search(r'M', cig) #find M in the cigar string. Split by soft?
#     return start_position


        #     if bioinfo.bitwise_strand_check(column)==FALSE:##strand is not reverse complement

        #     else: # if strand is reverse compelement 
        # else: 
        #     break 

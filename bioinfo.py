#!/usr/bin/env python
# Author: Lisa Adriani ladriani@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ACGTNagctn"
RNA_bases = "ACGUNagcun"

import re 

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return (ord(letter)) - 33

def qual_score(line: str):
    '''Converts the single character into a phred score and averages the phred score along the line.'''
    my_sum=0
    for letter in line:
        converted = convert_phred(letter)
        my_sum += converted 
    ave = my_sum / len(line)
    return ave

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # YOUR CODE HERE
    y = list()
    for x in range(101):
        y.append(0.0)
    return y

def validate_base_seq(DNA:str, RNA_flag: bool=False):
    '''This will validate if a string that is given is a base sequence, returning True or False'''
    DNA = DNA.upper()
    thelength= len(DNA)
    return thelength == DNA.count("A") + DNA.count("G") + DNA.count("C") + DNA.count("U" if RNA_flag else "T") +DNA.count("N")
    
def gc_content(DNAgc):
    '''This will return the proportion of GC content in a given DNA string'''
    DNAgc = DNAgc.upper()
    if validate_base_seq(DNAgc) == True:
        totGC = DNAgc.count("G") + DNAgc.count("C")
        GC = totGC/ len(DNAgc)
    return GC


def oneline_fasta(file, onelinefile):
    '''Takes a fasta file that has the sequence on multiple lines and rewrites the file with no new lines'''
    with open(onelinefile, "w") as o:
        with open(file, "r") as f:
            build_seq = ""
            i = 0
            for line in f:
                if line.startswith(">") == True:
                    header = line
                    if i > 0:
                        o.write(str(build_seq) + "\n")
                    o.write(str(header))
                    build_seq = ""
                if line.startswith(">") == False:
                    line = line.strip('\n') 
                    build_seq += line
                    # if i == 0:
                    #     o.write(str(build_seq) +"\n")
                    i += 1
            o.write(str(build_seq) + "\n")
    return

def revcomp(seq):
    '''Takes a string and returns the DNA reverse complement of it '''
    comp = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    bases = list(seq)
    rev = ""
    for base in bases:
        rev += comp[base]
    return rev[::-1]
        
def write_fastq_file(filename, array_name):
    '''takes a file name and a read number 1 or 2 and writes a fastq record (that is saved in an array) to the open file'''
    filename.write(array_name[0]+"\n"+array_name[1]+"\n"+array_name[2]+"\n"+array_name[3]+"\n")
    return

def bitwise_strand_check(POS1):
    '''Takes a bit flag and returns +/- where "-" is reverse complement, minus strand.'''
    bitW = int(POS1)
    if((bitW & 16)) == 16:
        position = "+"
    else:
        position = "-"
    return position

def get_UMI(POS0):
    '''Takes first column of SAM file and extracts the UMI'''
    splits=POS0.split(":")
    UMI = str(splits[7])
    return UMI

def correct_pos(cigar,strand,position):
    '''Takes the cigar string (as a string) and the position from the SAM file, and strand, and corrects for starting position'''
    cig = str(cigar)
    position = int(position)
    if strand == "+": # the plus strand
        findit = re.findall(r'(^\d+)([S]{1})',cig) #find the soft clipping at the beginning of the strand. i.e.[('20', 'S')]
        if findit != []: #if anything is found at all,
            soft = findit[0][0] # this will take the integer in front of the S -- data = [('20', 'S')]
            soft = int(soft)
            start_position = position - soft #subtract it from the position for the correct start position
        else: #if there is no left hand soft clipping
            start_position = position #position is your start position. 
    if strand == "-":
        findit = re.findall(r'(\d+)([A-Z]{1})',cig) # this will put data into tuples of # and letter i.e. 
        for letters in findit: 
            if letters[1] == "M":
                position = position + int(letters[0]) #adjust it!
            if letters[1] == "D":
                position = position + int(letters[0])
            if letters[1] == "N":
                position = position + int(letters[0])
        if findit[-1][1]== "S":
                position = position + int(findit[-1][0])
        start_position = position - 1 #subtract one because the position is inclusive. this is more technically correct. 

    return start_position


if __name__ == "__main__":
    assert convert_phred("E") == 36, "convert_phred working!"
    assert convert_phred("A") == 32, "convert_phred working!"

    assert qual_score("AEEA") == 34, "qual_score working!"
    assert qual_score("AEEEE") == 35.2, "qual_score working!"

    assert validate_base_seq(DNA_bases) == True, "validate_base_seq working!"
    assert validate_base_seq(RNA_bases, True) == True, "RNA validate_base_seq working!"
    assert validate_base_seq("IMSOTIRED", True) == False, "Validate_base_seq working"
    assert validate_base_seq("ATATATAT", True) == False, "RNA flag does work!"

    assert gc_content("GCGC") == 1, "gc_content on track!"
    assert gc_content("TATATA") == 0, "gc_content doing it, bestie!"
    assert gc_content("GATACG") == 0.5, "yeees!!"

    assert correct_left_soft("17S12M",20) == 3, "Soft Clip Position Adjustment working"
    #return
    # write tests for functions above
    # this if statement will ensure that the assert statements only run when you want to test it, and not in a jn for example.return
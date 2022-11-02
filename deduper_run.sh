#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=deduper.%j.out
#SBATCH --cpus-per-task=8
#SBATCH --job-name=deduper

conda activate bgmp_py310 

#/usr/bin/time -v samtools sort /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o sorted.sam

/usr/bin/time -v ./Adriani_deduper.py -f /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o test/bigone.sam -u unit_test_folder/STL96.txt 



##Testing Sam's Test files
#/usr/bin/time -v ./Adriani_deduper.py -f test.sam -o test/sams_output.sam -u unit_test_folder/STL96.txt 



#./Adriani_deduper.py --help
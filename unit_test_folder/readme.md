STL96.txt holds the 96 Known UMI sequences for our project. 

test1.sam has: 
    headers
    line 25/26--different cigar and position, with the soft clipping at the beginning of the line. Position should be corrected for. 
    line 27/28 -- all the same
    line 29/30 -- different cigar and position, soft clipping at the end of the line. Position should not be corrected for. 
    line 31/32 -- different strands
    line 33/34 -- different cigar strings, but the same position
    line 35/36 -- different UMIs, everything else the same
    line 37/38 -- same UMI and info as line 25 as an extra check at the end. 

output1.sam should have: 
    headers
    line 25, 27, 29, 31, 32, 33, 34, 35, 36. (I think)

    
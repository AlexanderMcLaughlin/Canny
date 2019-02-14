# Canny
A simple implementation of the Canny edge detection algorithm using convolution and a threshold value to determine where edges lie.

Note this is a Code Blocks project, however it can be compiled like so:

>gcc infile.pgm outfile1.pgm outfile2.pgm outfile3.pgm percentage

infile will be the input image, outfile1 will be the image's gradient magnitude, outfile2 will perform the peaks (topological mapping) algorithm on the image and show where the highest points are contiguous to other high points, outfile3 will find the outline of the image based on the percentage input, the more it is increased the less scrutinous the algorithm will be, should be between 60 and 90 for best results.

the input and output image files should be 256 pixels by 256 pixels for the program to do the entire image (this can be changed by updating the PICSIZE variable in canny.c.

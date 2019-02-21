#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Definitions of picture sizes
#define  PICSIZE 256
#define  MAXMASK 100

//Global picture variables
int pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];

double outpic1[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
double outpic3[PICSIZE][PICSIZE];

int    edgeflag[PICSIZE][PICSIZE];
double mask[MAXMASK][MAXMASK];

double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];

int histogram[PICSIZE];

double conv[PICSIZE][PICSIZE];

int main(int argc, char** argv)
{
    int     i, j, p, q, mr, centx, centy, areaOfTops=0, HI, LO, flag=0;
    double  maskval, sum1, sum2, sig, maxival, slope=0, percent=0;
    FILE    *fo1, *fo2, *fo3, *fp1;
    char    *foobar;

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get the input file
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get the first output file name
    foobar = *argv;
    fo1=fopen(foobar,"wb");

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get the second output file name
    foobar = *argv;
    fo2=fopen(foobar,"wb");

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get the third output file name
    foobar = *argv;
    fo3=fopen(foobar,"wb");

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get sigma value
    foobar = *argv;
    sig = atof(foobar);

    //Increment to the next command line argument
    argc--;
    argv++;

    //Get the percentage for part 4 (for threshold)
    foobar = *argv;
    percent = atof(foobar);


    //get the masked radius as sigma * 3
    mr = (int)(sig * 3);

    //Get the maximum mask value as 1/2 of the image pixel size
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    //Reading the header in just to move the file pointer to the correct position
    fscanf(fp1, "%*s\n%*s %*s\n%*s\n");

    //Read from the file into an array
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
          pic[i][j]= getc(fp1);
        }
    }

    //Get the values using the first derivative of the gaussian function
    for (p=-mr;p<=mr;p++)
    {
        for (q=-mr;q<=mr;q++)
        {
            //Get the x mask
            maskval = q*exp(-1*(((q*q)+(p*p))/(2*(sig*sig))));
            maskx[p+centy][q+centx] = maskval;

            //Get the y mask
            maskval = p*exp(-1*(((p*p)+(q*q))/(2*(sig*sig))));
            masky[p+centy][q+centx] = maskval;
        }
    }

    //Convolution with the gaussian function
    for (i=mr;i<=255-mr;i++)
    {
        for (j=mr;j<=255-mr;j++)
        {
            sum1= 0;
            sum2= 0;
            for (p=-mr;p<=mr;p++)
            {
                for (q=-mr;q<=mr;q++)
                {
                   sum1 += pic[i+p][j+q] * maskx[p+centy][q+centx];
                   sum2 += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    //Getting the gradient magnitude of the input image
    maxival = 0;
    for (i=mr;i<256-mr;i++)
    {
        for (j=mr;j<256-mr;j++)
        {
            outpic1[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) + (outpicy[i][j]*outpicy[i][j])));
            if (outpic1[i][j] > maxival)
                maxival = outpic1[i][j];
        }
    }

    //Output the header of a pgm format to the file
    fprintf(fo1,"P5\n256 256\n255\n");

    //Doing the scaling of the image afterwards and writing it to the file so there is no value overflow
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            outpic1[i][j] = (outpic1[i][j] / maxival) * 255;
            fprintf(fo1,"%c",(char)((int)(outpic1[i][j])));
        }
    }

    //Canny peaks tracing algorithm (similar to a topological map)
    for(i=mr;i<256-mr;i++)
    {
        for(j=mr;j<256-mr;j++)
        {
            if((outpicx[i][j]) == 0.0)
            {
                outpicx[i][j] = .00001;
            }

            slope = ((double)outpicy[i][j])/outpicx[i][j];

            if( (slope <= .4142)&&(slope > -.4142))
            {
                if((outpic1[i][j] > outpic1[i][j-1])&&(outpic1[i][j] > outpic1[i][j+1]))
                {
                    outpic2[i][j] = 255;
                }
            }
            else if( (slope <= 2.4142)&&(slope > .4142))
            {
                if((outpic1[i][j] > outpic1[i-1][j-1])&&(outpic1[i][j] > outpic1[i+1][j+1]))
                {
                    outpic2[i][j] = 255;
                }
            }
            else if( (slope <= -.4142)&&(slope > -2.4142))
            {
                if((outpic1[i][j] > outpic1[i+1][j-1])&&(outpic1[i][j] > outpic1[i-1][j+1]))
                {
                    outpic2[i][j] = 255;
                }
            }
            else
            {
                if((outpic1[i][j] > outpic1[i-1][j])&&(outpic1[i][j] > outpic1[i+1][j]))
                {
                    outpic2[i][j] = 255;
                }
            }
        }
    }

    //Output the pgm header to the file
    fprintf(fo2,"P5\n256 256\n255\n");

    //Writing it to the second file
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            fprintf(fo2,"%c",(char)((int)(outpic2[i][j])));
        }
    }


    //Algorithm to get the histogram of scaled magnitudes
    for(i=0; i<256; i++)
    {
        for(j=0; j<256; j++)
        {
            histogram[(int)outpic1[i][j]]++;
        }
    }

    //Get HI automatically
    for (HI=255; areaOfTops<((percent/100)*(256*256)) && HI>0; HI--)
        areaOfTops += histogram[HI];

    LO= 0.35*HI;

    //Output of the threshold
    printf("PERCENT=%.2lf\nHIGH = %d\nLOW = %d\n", percent,HI,LO);

    //Double thresholding
    for(i=0; i<256; i++)
    {
        for(j=0; j<256; j++)
        {
            if(outpic2[i][j]==255)
            {
                if(outpic1[i][j]>HI)
                {
                    outpic2[i][j]=0;
                    outpic3[i][j]=255;
                }
                else if(outpic1[i][j]<LO)
                {
                    outpic2[i][j]=0;
                    outpic3[i][j]=0;
                }

            }

        }

    }

    //Associates values below the threshold percentage but above the low
    //threshold value with an edge if the value is connected to a high value
    flag=1;
    while(flag)
    {
        flag=0;

        for(i=1; i<255; i++)
        {
            for(j=1; j<255; j++)
            {
                if(outpic2[i][j]==255)
                {
                    for(p=i-1; p<=i+1; p++)
                    {
                        for(q=j-1; q<=j+1; q++)
                        {
                            if(outpic3[i+p][j+q]==255)
                            {
                                outpic2[i][j]=0;
                                outpic3[i][j]=255;
                                flag=1;
                            }
                        }
                    }
                }
            }
        }
    }

    //Output the pgm header to the third file
    fprintf(fo3,"P5\n256 256\n255\n");

    //Writing it to the third file
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            fprintf(fo3,"%c",(char)((int)(outpic3[i][j])));
        }
    }

    return 0;
}

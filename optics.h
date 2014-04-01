#ifndef __OPTICS__
#define __OPTICS__

/********************************************************/
/*  Calculate optical ouptut using Jones matrix method  */
/*  and write to BMP                                    */
/********************************************************/

#include <ctime>
#include <time.h>
#include <float.h> 
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

//include LcSimParameters
#include "LcSimParameters.h"

/**********************************************************************/
/** Data structure/ function for writing bmp image                   **/
/**********************************************************************/
struct BMPHeader
{
    char bfType[2];       /* "BM" */
    int bfSize;           /* Size of file in bytes */
    int bfReserved;       /* set to 0 */
    int bfOffBits;        /* Byte offset to actual bitmap data (= 54) */
    int biSize;           /* Size of BITMAPINFOHEADER, in bytes (= 40) */
    int biWidth;          /* Width of image, in pixels */
    int biHeight;         /* Height of images, in pixels */
    short biPlanes;       /* Number of planes in target device (set to 1) */
    short biBitCount;     /* Bits per pixel (24 in this case) */
    int biCompression;    /* Type of compression (0 if no compression) */
    int biSizeImage;      /* Image size, in bytes (0 if no compression) */
    int biXPelsPerMeter;  /* Resolution in pixels/meter of display device */
    int biYPelsPerMeter;  /* Resolution in pixels/meter of display device */
    int biClrUsed;        /* Number of colors in the color table (if 0, use 
                             maximum allowed by biBitCount) */
    int biClrImportant;   /* Number of important colors.  If 0, all colors 
                             are important */
};

// Function to write BMP image
int write_bmp(const char *filename, int width, int height, char *rgb);



/**********************************************************************/
/** Data structure/ function for optics calcuation                   **/
/**********************************************************************/
struct optics_struct{
	double R,G,B;
};

typedef std::complex<double> dcomp;     //define complex number type

void matrix_mult(dcomp (&A)[2][2],dcomp (&B)[2][2],dcomp (&C)[2][2]);
void matrix_copy(dcomp (&A)[2][2],dcomp (&B)[2][2]);
void get_intensity(float *director,optics_struct (&intensity)[ISIZE*JSIZE]);
void BMP_write(int n,optics_struct (&intensity)[ISIZE*JSIZE]);
double dot(double (&A)[3],double (&B)[3]);
void cross(double (&A)[3],double (&B)[3], double (&C)[3]);


/**********************************************************************/
/** Mulitpy A*B and get C                                            **/
/**********************************************************************/
void matrix_mult(dcomp (&A)[2][2],dcomp (&B)[2][2],dcomp (&C)[2][2])
{
	dcomp val;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			val = dcomp(0.,0.);
			for (int k=0;k<2;k++)
			{
				val+=A[i][k]*B[k][j];
			}
			C[i][j]=val;
			
		}
	}
}// matrix_mult


/**********************************************************************/
/** Copy matrix A to matrix B                                        **/
/**********************************************************************/
void matrix_copy(dcomp (&A)[2][2],dcomp (&B)[2][2])
{
	for (int i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
		{
			A[i][j]=B[i][j];
		}
	}
}//end matrix copy


/**********************************************************************/
/** write bmp file                                                   **/
/**********************************************************************/
void BMP_write(int n,optics_struct (&intensity)[ISIZE*JSIZE])
{
	int k=0;
	char A[3*ISIZE*JSIZE];
	for (int i=0;i<ISIZE;i++)
	{
		for (int j=0;j<JSIZE;j++)
		{
			A[k]=char(255*intensity[i+ISIZE*j].R);
			k++;
			A[k]=char(255*intensity[i+ISIZE*j].G);
			k++;
			A[k]=char(255*intensity[i+ISIZE*j].B);
			k++;
		}
	}
	
	char filename[30];
	sprintf(filename,"BMP//PIC%d.bmp",n);
	write_bmp(filename,ISIZE,JSIZE,A);
}//BMP_write


/**********************************************************************/
/** calculate the intensity for the RGB channesl                     **/
/**********************************************************************/
void get_intensity(float *director,optics_struct (&intensity)[ISIZE*JSIZE])
{
	dcomp img,expipsiR,ExR,EyR,expipsiG,ExG,EyG,expipsiB,ExB,EyB;
	dcomp JoldR[2][2],JoldG[2][2],JoldB[2][2];
	dcomp JkR[2][2],JkG[2][2],JkB[2][2];
	dcomp JnewR[2][2],JnewG[2][2],JnewB[2][2];
	dcomp Eout1R[2][2],Eout1G[2][2],Eout1B[2][2];
	dcomp Polin[2][2]={{dcomp(1.,0.),dcomp(0.,0.)},{dcomp(0.,0.),dcomp(0.,0.)}};
	dcomp Polout[2][2]={{dcomp(0.,0.),dcomp(0.,0.)},{dcomp(0.,0.),dcomp(1.,0.)}};
	dcomp Ein[2][2] = {{dcomp(1.,0.),dcomp(0.,0.)},{dcomp(0.,0.),dcomp(1.,0.)}};
	img = dcomp (0., 1.); 

	double sintheta,costheta,deltan,neff,cosphi,sinphi,Nx,Ny,Nz,theta,phi,psiR,psiG,psiB;
	double sintheta2,costheta2,sincostheta,IntensityR,IntensityG,IntensityB;
	double nemax=NeMax,no=No;
	double lambdaR=0.650,lambdaG=0.510,lambdaB=0.475;
	double neno=nemax*no,ne2=nemax*nemax,no2=no*no;
	double dee = OPTICS_THICKNESS/double(KSIZE);
	
	for(int i=0;i<ISIZE;i++){
		for(int j=0;j<JSIZE;j++){
			//begin calc for one pixel
			   matrix_copy(JnewR,Ein);
			   matrix_copy(JnewG,Ein);
			   matrix_copy(JnewB,Ein);

			   for(int k=0;k<KSIZE;k++){
				   matrix_copy(JoldR,JnewR);
				   matrix_copy(JoldG,JnewG);
				   matrix_copy(JoldB,JnewB);

				   Nx = director[0+3*(i+ISIZE*(j+JSIZE*k))];
				   Ny = director[1+3*(i+ISIZE*(j+JSIZE*k))];
				   Nz = director[2+3*(i+ISIZE*(j+JSIZE*k))];

                   //handle Nz = 0 and Nx = zero edge cases
                   if(Nz*Nz<0.00001){Nz=0.001;}
                   if(Nx*Nx<0.00001){Nx=0.001;}
				   theta = atan(Ny/Nx);
		           phi = atan(sqrt(Nx*Nx+Ny*Ny)/Nz);
				        
		           sintheta=sin(theta);
	    	       costheta=cos(theta);
			       sintheta2=sintheta*sintheta;
				   costheta2=costheta*costheta;
				   sincostheta=sintheta*costheta;
				        
				   sinphi=sin(phi);
				   cosphi=cos(phi);
				        
				   neff= neno/sqrt(ne2*cosphi*cosphi+no2*sinphi*sinphi);
				   deltan=neff-no;
				   psiR=2.*PI*deltan*dee/lambdaR;
				   psiG=2.*PI*deltan*dee/lambdaG;
				   psiB=2.*PI*deltan*dee/lambdaB;
				        
				   expipsiR=exp(-img*psiR);
				   expipsiG=exp(-img*psiG);
				   expipsiB=exp(-img*psiB);

				   JkR[0][0]=costheta2+sintheta2*expipsiR;
				   JkR[0][1]=sincostheta*(1.-expipsiR);
				   JkR[1][0]=sincostheta*(1.-expipsiR);
				   JkR[1][1]=sintheta2+costheta2*expipsiR;

				   JkG[0][0]=costheta2+sintheta2*expipsiG;
				   JkG[0][1]=sincostheta*(1.-expipsiG);
				   JkG[1][0]=sincostheta*(1.-expipsiG);
				   JkG[1][1]=sintheta2+costheta2*expipsiG;

				   JkB[0][0]=costheta2+sintheta2*expipsiB;
				   JkB[0][1]=sincostheta*(1.-expipsiB);
				   JkB[1][0]=sincostheta*(1.-expipsiB);
				   JkB[1][1]=sintheta2+costheta2*expipsiB;

				   matrix_mult(JoldR,JkR,JnewR);
				   matrix_mult(JoldG,JkG,JnewG);
				   matrix_mult(JoldB,JkB,JnewB);
			   }

			   matrix_mult(Polin,JnewR,Eout1R);
			   matrix_mult(Polin,JnewG,Eout1G);
			   matrix_mult(Polin,JnewB,Eout1B);

			   ExR=Eout1R[0][0];
			   EyR=Eout1R[0][1];
			   ExG=Eout1G[0][0];
			   EyG=Eout1G[0][1];
			   ExB=Eout1B[0][0];
			   EyB=Eout1B[0][1];

			   IntensityR=real(EyR*conj(EyR));
			   IntensityG=real(EyG*conj(EyG));
			   IntensityB=real(EyB*conj(EyB));

			 //  printf("R = %f, G = %f B= %f\n",IntensityR,IntensityG,IntensityB);

			   intensity[i+ISIZE*j].R = IntensityR;
			   intensity[i+ISIZE*j].G = IntensityG;
			   intensity[i+ISIZE*j].B = IntensityB;
		}
	}

}//get_intensity


/*************************************************************************
    This function is from bmp.c

    Copyright (C) 2002 Hari Nair <hari@alumni.caltech.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
**************************************************************************/

int 
write_bmp(const char *filename, int width, int height, char *rgb)
{
    int i, j, ipos;
    int bytesPerLine;
    unsigned char *line;

    FILE *file;
    struct BMPHeader bmph;

    /* The length of each line must be a multiple of 4 bytes */

    bytesPerLine = (3 * (width + 1) / 4) * 4;

    strcpy(bmph.bfType, "BM");
    bmph.bfOffBits = 54;
    bmph.bfSize = bmph.bfOffBits + bytesPerLine * height;
    bmph.bfReserved = 0;
    bmph.biSize = 40;
    bmph.biWidth = width;
    bmph.biHeight = height;
    bmph.biPlanes = 1;
    bmph.biBitCount = 24;
    bmph.biCompression = 0;
    bmph.biSizeImage = bytesPerLine * height;
    bmph.biXPelsPerMeter = 0;
    bmph.biYPelsPerMeter = 0;
    bmph.biClrUsed = 0;       
    bmph.biClrImportant = 0; 

    file = fopen (filename, "wb");
    if (file == NULL) return(0);
  
    fwrite(&bmph.bfType, 2, 1, file);
    fwrite(&bmph.bfSize, 4, 1, file);
    fwrite(&bmph.bfReserved, 4, 1, file);
    fwrite(&bmph.bfOffBits, 4, 1, file);
    fwrite(&bmph.biSize, 4, 1, file);
    fwrite(&bmph.biWidth, 4, 1, file);
    fwrite(&bmph.biHeight, 4, 1, file);
    fwrite(&bmph.biPlanes, 2, 1, file);
    fwrite(&bmph.biBitCount, 2, 1, file);
    fwrite(&bmph.biCompression, 4, 1, file);
    fwrite(&bmph.biSizeImage, 4, 1, file);
    fwrite(&bmph.biXPelsPerMeter, 4, 1, file);
    fwrite(&bmph.biYPelsPerMeter, 4, 1, file);
    fwrite(&bmph.biClrUsed, 4, 1, file);
    fwrite(&bmph.biClrImportant, 4, 1, file);
  
    line =(unsigned char*) malloc(bytesPerLine);
    if (line == NULL)
    {
        fprintf(stderr, "Can't allocate memory for BMP file.\n");
        return(0);
    }

    for (i = height - 1; i >= 0; i--)
    {
        for (j = 0; j < width; j++)
        {
            ipos = 3 * (width * i + j);
            line[3*j] = rgb[ipos + 2];
            line[3*j+1] = rgb[ipos + 1];
            line[3*j+2] = rgb[ipos];
        }
        fwrite(line, bytesPerLine, 1, file);
    }

    free(line);
    fclose(file);

    return(1);
}//writebmp




#endif
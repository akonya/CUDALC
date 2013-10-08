/**********************************************************************/
/** Simulation of cholestaric with variable q                        **/
/** Version 1.0                                                      **/
/** January 10, 2012                                                 **/
/** Written by: Andrew Konya                                         **/
/**                                                                  **/
/**********************************************************************/

#include <ctime>
#include <time.h>
#include <float.h> 
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>




/**********************************************************************/
/** Define constants to be used in the simulation                    **/
/**********************************************************************/
#define PI              3.141592654        //the constant

#define Nsteps          100000             //number of times to itterate the system
#define dt              0.04               //size of each step
#define Xsize           150                //size of system in X direction
#define Ysize           45                 //size of system in Y direction
#define Zsize           10                 //size of system in Z direction
#define NP              Xsize*Ysize*Zsize  //number of voxels in the system

#define NeMax           1.7                //max refractive index
#define No              1.5                //ordinary refractive index
#define p               15.0               //charachteristic pitch  (I think in units of voxels)
#define q0              2.0*PI/p           //max q of system
#define d               3.5                //cell thicknes in microns


/**********************************************************************/
/** constatns/functions for random number generator                  **/
/** DONT CHANGE THESE                                                **/
/**********************************************************************/
#define MT_LEN          624
#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

int mt_index;
unsigned long mt_buffer[MT_LEN];

void mt_init();
unsigned long randmt();
double genrand();
void purge();
double randgauss(double sigma, double mean);



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

int write_bmp(const char *filename, int width, int height, char *rgb);





/**********************************************************************/
/** Define data structures to be used in system                      **/
/**                                                                  **/
/** Examples: --To get the torque in the x direction at point i,j,k  **/
/**             you would type: vox[i][j][k].tx                      **/
/**                                                                  **/
/**			  -- To get the q at point i,j, you would type:          **/
/**                           vox[i][j][k].q                         **/
/**********************************************************************/
struct vox_struct{
	double nx,ny,nz,tx,ty,tz,E;
};

vox_struct vox[Xsize][Ysize][Zsize];


struct optics_struct{
	double R,G,B;
};

optics_struct intensity[Xsize][Ysize];

typedef std::complex<double> dcomp;     //define complex number type









/**********************************************************************/
/** Declare functions to be used, the bodies of the functions are    **/
/** below the main loop                                              **/
/**********************************************************************/
void initialize(vox_struct (&vox)[Xsize][Ysize][Zsize]);
float calc_torque(vox_struct (&vox)[Xsize][Ysize][Zsize],double KbT);
void zero_torque(vox_struct (&vox)[Xsize][Ysize][Zsize]);
void update_director(vox_struct (&vox)[Xsize][Ysize][Zsize]);
void matrix_mult(dcomp (&A)[2][2],dcomp (&B)[2][2],dcomp (&C)[2][2]);
void matrix_copy(dcomp (&A)[2][2],dcomp (&B)[2][2]);
void get_intensity(vox_struct (&vox)[Xsize][Ysize][Zsize],optics_struct (&intensity)[Xsize][Ysize]);
void BMP_write(int n,optics_struct (&intensity)[Xsize][Ysize]);
void printVTKframe(vox_struct (&vox)[Xsize][Ysize][Zsize],int step);



/**********************************************************************/
/** Main loop of the system                                          **/
/**********************************************************************/
int main()
{
	srand(98237);    //seed random number generator
	mt_init();       //initialize random number generator
	purge();         //free up memory in random number generator

	double KbT=0.0;
	double Etot;
  int switcher = 0;

	initialize(vox);
  printVTKframe(vox,0);

	for (int step=1;step<Nsteps;step++){
    
    //temperature control
    if(KbT>0.05 && switcher == 0){
      KbT = KbT*0.9997;
    }else{
      switcher = 1;
    //  KbT += 0.000001;
    }//if KbT

		Etot = calc_torque(vox,KbT);          
		update_director(vox);          
		zero_torque(vox);

		if (step%100==0){     //we can check things every 100 steps
			printf("%d of %d KbT = %f Energy = %f\n",step,Nsteps,KbT,Etot);
			printVTKframe(vox,step);
      get_intensity(vox,intensity);
			BMP_write(step,intensity);
		}
	}
	return 0;
}





/**********************************************************************/
/** Function bodies                                                  **/
/**********************************************************************/

/**********************************************************************/
/** write directors to VTK frame                                     **/
/**********************************************************************/
void printVTKframe(vox_struct (&vox)[Xsize][Ysize][Zsize],int step){
  int npoints;
  int count;
  double nnx,nny,nnz;

  //count total points
  count = 0;
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        count++;    
      }//fork
    }//forj
  }//fori
  npoints=count;
 
  //write filename and open file to write to
  char fout[60];
  sprintf(fout,"VTK//dir%d.vtk",step);
  FILE*out;
  out = fopen(fout,"w");
  
  //write VTK header
  fprintf(out,"# vtk DataFile Version 3.1\n");
  fprintf(out,"director profile\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(out,"\n");
  fprintf(out,"POINTS %d FLOAT\n",npoints);
  
  //loop over all points to get postions
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        fprintf(out,"%f %f %f\n",float(i),float(j),float(k));      
      }//fork
    }//forj
  }//fori
  fprintf(out,"\n");

  //cells
  fprintf(out,"CELLS %d %d\n",npoints,npoints*2);
  count = 0;
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        fprintf(out,"1 %d\n",count);
        count++;    
      }//fork
    }//forj
  }//fori
  fprintf(out,"\n");

  //cell types
  fprintf(out,"CELL_TYPES %d\n",npoints);
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        fprintf(out,"1\n");
      }//fork
    }//forj
  }//fori
  fprintf(out,"\n");

  //director data
  fprintf(out,"POINT_DATA %d\n",npoints);
  fprintf(out,"VECTORS director FLOAT\n");
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        nnx = vox[i][j][k].nx;
        nny = vox[i][j][k].ny;
        nnz = vox[i][j][k].nz;
        fprintf(out,"%f %f %f\n",nnx,nny,nnz);
      }//fork
    }//forj
  }//fori
  fprintf(out,"\n");

  //ENERGY
  fprintf(out,"SCALARS energy FLOAT 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
        fprintf(out,"%f\n",vox[i][j][k].E);
      }//fork
    }//forj
  }//fori
  fprintf(out,"\n");


  fclose(out);
}//printVTKframe




/**********************************************************************/
/** Make a picture                                                   **/
/**********************************************************************/
void BMP_write(int n,optics_struct (&intensity)[Xsize][Ysize])
{
	int k=0;
	char A[3*Xsize*Ysize];
	for (int i=0;i<Xsize;i++)
	{
		for (int j=0;j<Ysize;j++)
		{
			A[k]=char(255*intensity[i][j].R);
			k++;
			A[k]=char(255*intensity[i][j].G);
			k++;
			A[k]=char(255*intensity[i][j].B);
			k++;
		}
	}
	
	
	
	char filename[30];
	sprintf(filename,"pics//PIC%d.bmp",n);
	write_bmp(filename,Ysize,Xsize,A);
}


/**********************************************************************/
/** calculate the intensity for the RGB channesl                     **/
/**********************************************************************/
void get_intensity(vox_struct (&vox)[Xsize][Ysize][Zsize],optics_struct (&intensity)[Xsize][Ysize])
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
	double dee = d/double(Zsize);
	
	for(int i=0;i<Xsize;i++){
		for(int j=0;j<Ysize;j++){
			//begin calc for one pixel
			   matrix_copy(JnewR,Ein);
			   matrix_copy(JnewG,Ein);
			   matrix_copy(JnewB,Ein);

			   for(int k=0;k<Zsize;k++){
				   matrix_copy(JoldR,JnewR);
				   matrix_copy(JoldG,JnewG);
				   matrix_copy(JoldB,JnewB);

				   Nx=vox[i][j][k].nx;
				   Ny=vox[i][j][k].ny;
				   Nz=vox[i][j][k].nz;
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

			   intensity[i][j].R = IntensityR;
			   intensity[i][j].G = IntensityG;
			   intensity[i][j].B = IntensityB;
		}
	}

}

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
}


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
/** Update the director                                              **/
/**********************************************************************/
void update_director(vox_struct (&vox)[Xsize][Ysize][Zsize])
{
	double wx,wy,wz,nnx,nny,nnz,nix,niy,niz,norm;
	double f = 1.0;

	for (int i=1;i<Xsize-1;i++){             //i is position in x direction
		for(int j=1;j<Ysize-1;j++){          //j is position in y direction
			for(int k=1;k<Zsize-1;k++){      //k is position in z direction

				nix = vox[i][j][k].nx;
				niy = vox[i][j][k].ny;
				niz = vox[i][j][k].nz;

				wx = f*vox[i][j][k].tx;
				wy = f*vox[i][j][k].ty;
				wz = f*vox[i][j][k].tz;

				nnx = nix+(wy*niz-wz*niy)*dt;
				nny = niy+(wz*nix-wx*niz)*dt;
				nnz = niz+(wx*niy-wy*nix)*dt;

				//printf("nnx = %f nny = %f nnz=%f\n",nnx,nny,nnz);

				norm = sqrt(nnx*nnx+nny*nny+nnz*nnz);
				vox[i][j][k].nx = nnx/norm;
				vox[i][j][k].ny = nny/norm;
				vox[i][j][k].nz = nnz/norm;
				
			}
		}
	}
}




/**********************************************************************/
/** Zero out the toruqe at each voxel                                **/
/**********************************************************************/
void zero_torque(vox_struct (&vox)[Xsize][Ysize][Zsize])
{
	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction
				vox[i][j][k].tx = 0.0;
				vox[i][j][k].ty = 0.0;
				vox[i][j][k].tz = 0.0;
			}
		}
	}
}



/**********************************************************************/
/** Calculate the toqrue on each director and store it               **/
/** no thermostate in yet                                            **/
/**********************************************************************/
float calc_torque(vox_struct (&vox)[Xsize][Ysize][Zsize],double KbT)
{
	int ip1,jp1,kp1;
	double dott,dott2,crossx,crossy,crossz,nix,niy,niz,njx,njy,njz,Tx,Ty,Tz,qSUM,xkick,ykick,zkick;
	double sig = sqrt(41.2*KbT);
	double Etot=0.0,Eloc;

	for (int i=0;i<Xsize;i++){							//i is position in x direction
		if (i+1==Xsize){ip1 = 0;}else{ip1 = i+1;}       //periodic boundary conditions in x
		for(int j=0;j<Ysize;j++){					    //j is position in y direction
			if (j+1==Ysize){jp1 = 0;}else{jp1 = j+1;}   //periodic boundary conditions in y
			for(int k=0;k<Zsize-1;k++){					//k is position in z direction
				kp1 = k+1;

        Eloc = 0.0;

				nix = vox[i][j][k].nx;
				niy = vox[i][j][k].ny;
				niz = vox[i][j][k].nz;

				//Neighbor in +x (ip1) direction

				xkick = randgauss(sig,0.0);
				ykick = randgauss(sig,0.0);
				zkick = randgauss(sig,0.0);

				njx = vox[ip1][j][k].nx;
				njy = vox[ip1][j][k].ny;
				njz = vox[ip1][j][k].nz;

				dott = nix*njx+niy*njy+niz*njz;
				dott2 = dott*dott;
				crossx = niy*njz-niz*njy;
				crossy = niz*njx-nix*njz;
				crossz = nix*njy-niy*njx;

				Tx = crossx*dott-q0*(crossx*crossx+dott*nix*njx-dott2)+xkick;
				Ty = crossy*dott-q0*(crossx*crossy+dott*nix*njy)+ykick;
				Tz = crossz*dott-q0*(crossx*crossz+dott*nix*njz)+zkick;

				vox[i][j][k].tx +=Tx;
				vox[i][j][k].ty +=Ty;
				vox[i][j][k].tz +=Tz;

				vox[ip1][j][k].tx +=-Tx;
				vox[ip1][j][k].ty +=-Ty;
				vox[ip1][j][k].tz +=-Tz;

				Eloc += 1.0+ q0*(dott*crossx)- dott*dott;        //calculate total elastic energy
        
				//Neighbor in +y (jp1) direction

				xkick = randgauss(sig,0.0);
				ykick = randgauss(sig,0.0);
				zkick = randgauss(sig,0.0);

				njx = vox[i][jp1][k].nx;
				njy = vox[i][jp1][k].ny;
				njz = vox[i][jp1][k].nz;

				dott = nix*njx+niy*njy+niz*njz;
				dott2 = dott*dott;
				crossx = niy*njz-niz*njy;
				crossy = niz*njx-nix*njz;
				crossz = nix*njy-niy*njx;

				Tx = crossx*dott-q0*(crossx*crossy+dott*niy*njx)+xkick;
				Ty = crossy*dott-q0*(crossy*crossy+dott*niy*njy-dott2)+ykick;
				Tz = crossz*dott-q0*(crossz*crossy+dott*niy*njz)+zkick;

				vox[i][j][k].tx +=Tx;
				vox[i][j][k].ty +=Ty;
				vox[i][j][k].tz +=Tz;

				vox[i][jp1][k].tx +=-Tx;
				vox[i][jp1][k].ty +=-Ty;
				vox[i][jp1][k].tz +=-Tz;

				Eloc += 1.0+ q0*(dott*crossy)- dott*dott;        //calculate total elastic energy

				//Neighbor in +z (kp1) direction

				xkick = randgauss(sig,0.0);
				ykick = randgauss(sig,0.0);
				zkick = randgauss(sig,0.0);

				njx = vox[i][j][kp1].nx;
				njy = vox[i][j][kp1].ny;
				njz = vox[i][j][kp1].nz;

				dott = nix*njx+niy*njy+niz*njz;
				dott2 = dott*dott;
				crossx = niy*njz-niz*njy;
				crossy = niz*njx-nix*njz;
				crossz = nix*njy-niy*njx;

				Tx = crossx*dott-q0*(crossx*crossz+dott*niz*njx)+xkick;
				Ty = crossy*dott-q0*(crossy*crossz+dott*niz*njy)+ykick;
				Tz = crossz*dott-q0*(crossz*crossz+dott*niz*njz-dott2)+zkick;

				vox[i][j][k].tx +=Tx;
				vox[i][j][k].ty +=Ty;
				vox[i][j][k].tz +=Tz;

				vox[i][j][kp1].tx +=-Tx;
				vox[i][j][kp1].ty +=-Ty;
				vox[i][j][kp1].tz +=-Tz;

				Eloc += 1.0+ q0*(dott*crossz)- dott*dott;        //calculate total elastic energy
        Etot += Eloc;
        vox[i][j][k].E = Eloc;
			}
		}
	}
	return (Etot/float(NP));    //function returns total elastic energy
}



/**********************************************************************/
/** Sets the initial director configuration of the system as well    **/
/** as the initial q's                                               **/
/**   --for now q's will be +-0.5 for bulk and 0.0 of substrates     **/
/**********************************************************************/
void initialize(vox_struct (&vox)[Xsize][Ysize][Zsize])
{
	double theta,phi,nx,ny,nz,norm;

	for (int i=0;i<Xsize;i++){             //i is position in x direction
		for(int j=0;j<Ysize;j++){          //j is position in y direction
			for(int k=0;k<Zsize;k++){      //k is position in z direction

        nx = 0.0;
        ny = 0.0;
        nz = 0.0;


       if(i==0 || i==(Xsize-1)){
            nx = 1.0;
            ny = 0.0001;
            nz = 0.0001;
        }//if i
       if(j==0 || j==(Ysize-1)){
            nx = 0.0001;
            ny = 1.0;
            nz = 0.0001;
        }//if j
        if(k==0 || k==(Zsize-1)){
            nx = 0.0001;
            ny = 0.0001;
            nz = 1.0;
        }//if i
        norm = sqrt(nx*nx+ny*ny+nz*nz);
        if(norm<0.1){
          nx = genrand()-0.5;//1.0/sqrt(2.0);
          ny = genrand()-0.5;//1.0/sqrt(2.0);
          nz = genrand()-0.5;//0.0;
          norm = sqrt(nx*nx+ny*ny+nz*nz);
          nx = nx/norm;
          ny = ny/norm;
          nz = nz/norm;

        }else{
          nx = nx/norm;
          ny = ny/norm;
          nz = nz/norm;
        }

        
				//store director directions
				vox[i][j][k].nx = nx;//sin(theta)*cos(phi);
				vox[i][j][k].ny = ny;//sin(theta)*sin(phi);
				vox[i][j][k].nz = nz;//cos(theta);

				//zero out torque
				vox[i][j][k].tx = 0.0;
				vox[i][j][k].ty = 0.0;
				vox[i][j][k].tz = 0.0;
				
			}
		}
	}
}




/**********************************************************************/
/** All function bodies for random number generator using MT         **/
/**********************************************************************/

void mt_init() {
    int i;
    for (i = 0; i < MT_LEN; i++)
        mt_buffer[i] = rand();
    mt_index = 0;
}



unsigned long randmt() {
    unsigned long * b = mt_buffer;
    int idx = mt_index;
    unsigned long s;
    int i;
	
    if (idx == MT_LEN*sizeof(unsigned long))
    {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < MT_LEN-1; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }
        
        s = TWIST(b, MT_LEN-1, 0);
        b[MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
    mt_index = idx + sizeof(unsigned long);
    return *(unsigned long *)((unsigned char *)b + idx);
}
double genrand()
{
unsigned long b=randmt();
(double)b;
double a;
a=b*0.000000000232830644;
return a;
}
void purge()
{
for (int j = 0; j < 2500; j++)
	{
		genrand();
	}
}

/**********************************************************************/
/** Generates numbers from a normal distrobution                     **/
/**********************************************************************/

double randgauss(double sigma, double mean)
{
float x1, x2, w, y1, y2, gauss;
 
         do {
                 x1 = 2.0 * genrand() - 1.0;
                 x2 = 2.0 * genrand() - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
		gauss = y2*sigma + mean;
return gauss;
}


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
}

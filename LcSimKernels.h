#ifndef __LCSIMKERNELS__
#define __LCSIMKERNELS__

#include "main.h"

//device vector opperations

//dot product
__device__ float dot(float (&A)[3], float (&B)[3]){
  return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
}//dot

//crossproduct
__device__ void cross(float (&A)[3], float (&B)[3], float (&C)[3]){
  C[0] = A[1]*B[2]-A[2]*B[1];
  C[1] = A[2]*B[0]-A[0]*B[2];
  C[2] = A[0]*B[1]-A[1]*B[0];
}//cross



//Kernel to Calculate torque
__global__ void calculateTorqueKernel(float *director_d
                                    , float *torque_d
                                    , int iSize
                                    , int jSize
                                    , int kSize
                                    , float q0){
  //declerations
  int ib, jb, kb;
  float naDOTnb, naDOTnb2, naDOTrab, nbDOTrab, naCROSSnbDOTrab;
  float naCROSSnb[3];
  float na[3];
  float nb[3];
  float TqaSum[3] = {0.0};
  float Tqa,Tqb,Tqc;
  float rab[3];
  int di[3] = {1,0,0};
  int dj[3] = {0,1,0};
  int dk[3] = {0,0,1}; 
  int shift[3] = {1,2,3}; //shift for storing torques to be summed 
 
  //get thread id
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //compute 3d indicies from 1d index
  int ka = tid/(iSize*jSize);
  int ja = (tid-(ka*iSize*jSize))/iSize;
  int ia = tid - iSize*(ja+ka*jSize);
 
  //only continue thread if ...
  if(ia<iSize-1 && ja<jSize-1 && ka<kSize-1){ 

    //get na.. looping over x,y,z
    for(int cord=0;cord<3;cord++){
      na[cord] = director_d[cord+3*(ia+iSize*(ja+jSize*ka))];
    }//cord

    //loop over neighbor particles
    for(int b=0;b<3;b++){

      //get rab vector
      rab[0] = float(di[b]);
      rab[1] = float(dj[b]);
      rab[2] = float(dk[b]);

      //get indiceis of neighbor 
      ib = ia+di[b];
      jb = ja+dj[b];
      kb = ka+dk[b];

      //get nb.. looping over x,y,z
      for(int cord=0;cord<3;cord++){
        nb[cord] = director_d[cord+3*(ib+iSize*(jb+jSize*kb))];
      }//cord

     //terms for torque calculaiton
     naDOTnb = dot(na,nb);
     naDOTnb2 = naDOTnb*naDOTnb;
     naDOTrab = dot(na,rab);
     nbDOTrab = dot(nb,rab);
     cross(na,nb,naCROSSnb);
     naCROSSnbDOTrab = dot(naCROSSnb,rab);

     //calculate torque compoents 
     for(int cord=0;cord<3;cord++){
     
       //common toqrue component
       Tqc= 2.0*naDOTnb*naCROSSnb[cord] \
                   + q0*( naCROSSnbDOTrab*naCROSSnb[cord] \
                   + naDOTnb2*rab[cord]);

       //na unique trq component
       Tqa = q0*naDOTnb*naDOTrab*nb[cord];

       //nb unique trq component (minus sign from rab -> rba)
       Tqb = -q0*naDOTnb*nbDOTrab*na[cord];

       //add trq for a component
       TqaSum[cord] += Tqc - Tqa;

       //add trq compoment for nb to global memeory
       torque_d[shift[b]+4*(cord+3*(ia+iSize*(ja+jSize*ka)))] = -Tqc-Tqb;

     }//cord 

    }//for b

  }//if (tid< ... )

}//

#endif

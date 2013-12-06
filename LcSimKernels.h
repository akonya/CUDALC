#ifndef __LCSIMKERNELS__
#define __LCSIMKERNELS__

#include "main.h"

//Kernel to Calculate torque
__global__ void calculateTorqueKernel(float *director_d
                                    , float *torque_d
                                    , int iSize
                                    , int jSize
                                    , int kSize){
  //get thread id
  int tid = threadId.x + blockIdx.x * blockDim.x;
  int i,j,k;

  //only continue thread if live site
  if(tid<iSize*jSize*kSize){
    //compute 3D indicies from 1D index
    k = tid/(iSize*jSize);
    j = (tid-(k*iSize*jSize))/iSize;
    i = tid - iSize*(j+k*jSize);

  }//if (tid< ... )

}//

#endif

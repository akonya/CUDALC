#ifndef __LCSIM__
#define __LCSIM__

//standard libraries
#include <stdio.h>
#include <math.h>

//include Kernel functions
#include "LcSimKernels.h"
#include "LcSimParameters.h"

//include files to print from Kernel
#include "cuPrintf.cu"
#include "cuPrintf.cuh"

//include curand
#include <curand_kernel.h>

//  Class for LC simulation
class LcSim{
  //class varriables
  public:
    float *director;
    float *energy;
    float *director_d;
    float *energy_d;
    float *torque_d;
    int iSize_,jSize_,kSize_;
    int threadsPerBlock;
    int blocksPerKernel;
    curandState *states_d;

  //class constructor/destructor
  LcSim();
  ~LcSim();

  //director initialization method
  void initDirector();
  //GPU initialization method
  void initDevice();
  //copy host data to device (director)
  void sendDataToDevice();
  //get data from device
  void getDataFromDevice();
  //calculuate torque
  void calculateTorque(float KbT); //launches calculateTorqueKernel
  //update director
  void updateDirector(); //launches updateDirectorKernel
  //print VTK frame 
  void printVtkFrame(int step);
  //clean up memory after simulation complete
  void shutdown();
  //simple random number [0,1]
  float randf();

};

//LcSim constructor
LcSim::LcSim(){
  //set sized of cubic grid
  iSize_ = ISIZE;
  jSize_ = JSIZE;
  kSize_ = KSIZE;
  //initialize director and energy arrays
  director = new float[3*iSize_*jSize_*kSize_];
  energy = new float[iSize_*jSize_*kSize_];
}//LcSim

//LcSim deconstructor
LcSim::~LcSim(){
  delete director;
  director = NULL;
  delete energy;
  energy = NULL;  
}//LcSim

//initilize LC director field
void LcSim::initDirector(){
  //declare
  float nx, ny, nz, norm;

  //note: director indexed by the convention
  //      director[i][j][k][cord] -> director[cord+3*(i+iSize_*(j+jSize_*k))]
  for (int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        nx = 2.0*(randf()-0.5);
        ny = 2.0*(randf()-0.5);
        nz = 2.0*(randf()-0.5);
        norm = sqrt(nx*nx+ny*ny+nz*nz);
       
        director[0+3*(i+iSize_*(j+jSize_*k))] = nx/norm;
        director[1+3*(i+iSize_*(j+jSize_*k))] = ny/norm;
        director[2+3*(i+iSize_*(j+jSize_*k))] = nz/norm;
      }//k
    }//j
  }//i                     
}//initializeDirectori

//Initialize memory and parametrs for the GPU
void LcSim::initDevice(){
  //set threadsPerBlock (temoporarlily here to DEBUG)
  threadsPerBlock=128;

  //calculate blocks to excicute
  blocksPerKernel = (iSize_*jSize_*kSize_+threadsPerBlock)/threadsPerBlock;

  //allocate memory for director on GPU
  cudaMalloc( (void**) &director_d
            , 3*iSize_*jSize_*kSize_*sizeof(float));
  //allocate memory for energy on GPU
  cudaMalloc( (void**) &energy_d
            , iSize_*jSize_*kSize_*sizeof(float));
  //allocate memory for torque on GPU
  cudaMalloc( (void**) &torque_d
            , 4*3*iSize_*jSize_*kSize_*sizeof(float));
  //set inital torque to zeros
  cudaMemset(torque_d,0,4*3*iSize_*jSize_*kSize_*sizeof(float));

  //initilize cudaPrintf
  cudaPrintfInit();

  //intialize memeory for curandStates
  cudaMalloc((void**)&states_d, iSize_*jSize_*kSize_*sizeof(curandState));

  //run initializaiton kernel
  setupCurandKernel<<<blocksPerKernel,threadsPerBlock>>>(states_d,iSize_*jSize_*kSize_);

}//initializeGPU

//send data from host to device (director)
void LcSim::sendDataToDevice(){
  //copy director on host to device
  cudaMemcpy( director_d
            , director
            , 3*iSize_*jSize_*kSize_*sizeof(float)
            , cudaMemcpyHostToDevice);
}//sendDataToDevice

//get data from device - director and energy
void LcSim::getDataFromDevice(){
  //copy director on device to host
  cudaMemcpy( director
            , director_d
            , 3*iSize_*jSize_*kSize_*sizeof(float)
            , cudaMemcpyDeviceToHost );
  //copy energy on device to host
  cudaMemcpy( energy
            , energy_d
            , iSize_*jSize_*kSize_*sizeof(float)
            , cudaMemcpyDeviceToHost );
}//getDataFromDevice

//calculate toqrue using GPU (lanunch GPU torque kernel)
void LcSim::calculateTorque(float KbT){
  //launch caclulateTorqueKernel on GPU
  calculateTorqueKernel<<<blocksPerKernel,threadsPerBlock>>>(
             director_d
           , torque_d
           , iSize_
           , jSize_
           , kSize_ 
           , QZERO
           , sqrt(41.2*KbT)
           , states_d);
 
  //print buffer from cuPrintf
  cudaPrintfDisplay(stdout,true);
}//calculate torque

//update director positions
void LcSim::updateDirector(){
  //launch updateDirectorKernel on GPU
  updateDirectorKernel<<<blocksPerKernel,threadsPerBlock>>>(
             director_d
           , torque_d
           , iSize_
           , jSize_
           , kSize_ 
           , DELTAT);

  //print buffer from cuPrintf
  cudaPrintfDisplay(stdout,true);
}//updateDirector

//Print VTK frame into /VTK file
void LcSim::printVtkFrame(int step){

  printf("step = %d\n",step);  

  //local delclerations
  int count;
  float nx,ny,nz;
  
  //calculate total sites
  int totalPoints = iSize_*jSize_*kSize_;

  //write filename and open file to write
  char fout[60];
  sprintf(fout, "VTK//dir%d.vtk",step);
  FILE*out;
  out = fopen(fout,"w");

  //write VKT header
  fprintf(out,"# vtk DataFile Version 3.1\n");
  fprintf(out,"director profile\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(out,"\n");

  //write point locations
  fprintf(out," POINTS %d FLOAT\n",totalPoints);
  for(int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        fprintf(out,"%f %f %f\n",float(i),float(j),float(k)); 
      }//k
    }//j
  }//i
  fprintf(out,"\n");

  //write cell
  fprintf(out,"CELLS %d %d\n",totalPoints,2*totalPoints);
  count = 0;
  for(int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        fprintf(out,"1 %d\n",count);
        count++; 
      }//k
    }//j
  }//i
  fprintf(out,"\n");

  //write cell types
  fprintf(out,"CELL_TYPES %d\n",totalPoints);
  for(int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        fprintf(out,"1\n"); 
      }//k
    }//j
  }//i
  fprintf(out,"\n");

  //write director data
  fprintf(out,"POINT_DATA %d\n",totalPoints);
  fprintf(out,"VECTORS director FLOAT\n");
  for(int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        nx =  director[0+3*(i+iSize_*(j+jSize_*k))];         
        ny =  director[1+3*(i+iSize_*(j+jSize_*k))];         
        nz =  director[2+3*(i+iSize_*(j+jSize_*k))];         
        fprintf(out,"%f %f %f\n",nx,ny,nz); 
      }//k
    }//j
  }//i
  fprintf(out,"\n");

  //write ENERGY at each lattice site
  fprintf(out,"SCALARS energy FLOAT 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for(int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        // get energy by "energy[i+iSize_*(j+jSize_*k)]" once calcualted       
        fprintf(out,"%f\n",0.0); 
      }//k
    }//j
  }//i
  fprintf(out,"\n");

}//printVtkFrame

//clean up memeory on device
void LcSim::shutdown(){
  cudaFree(director_d);
  cudaFree(torque_d);
  cudaFree(states_d);
}//shutdowni

//random number (0,1)
float LcSim::randf(){
  return float(rand())/float(RAND_MAX);
}//randf
#endif

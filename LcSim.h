#ifndef __LCSIM__
#define __LCSIM__

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
  //class constructor/destructor
  LcSim(int iSize,int jSize, int kSize);
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
  void calculateTorque(); //launches calculateTorqueKernel

};

//LcSim constructor
LcSim::LcSim(int iSize, int jSize, int kSize){
  //set sized of cubic grid
  iSize_ = iSize;
  jSize_ = jSize;
  kSize_ = kSize;
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
  //note: director indexed by the convention
  //      director[i][j][k][cord] -> director[cord+3*(i+iSize_*(j+jSize_*k))]
  for (int i=0;i<iSize_;i++){
    for(int j=0;j<jSize_;j++){
      for(int k=0;k<kSize_;k++){
        director[0+3*(i+iSize_*(j+jSize_*k))] = 1.0;
        director[1+3*(i+iSize_*(j+jSize_*k))] = 0.0;
        director[2+3*(i+iSize_*(j+jSize_*k))] = 0.0;
      }//k
    }//j
  }//i                     
}//initializeDirectori

//Initialize memory and parametrs for the GPU
void LcSim::initDevice(){
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
void LcSim::calculateTorque(){
 /*
  calculateTorqueKernel<<<blocksPerKernel,threadsPerBlock>>>(
             director_d
           , torque_d);
*/
}//calculate torque
#endif

#ifndef __LCSIM__
#define __LCSIM__

//  Class for LC simulation
class LcSim{
  //class varriables
  public:
    float *director;
    float *energy;
    float *gpuDirector;
    float *gpuEnergy;
    float *gpuTorque;
    int iSize_,jSize_,kSize_;
  //class constructor/destructor
  LcSim(int iSize,int jSize, int kSize);
  ~LcSim();
  //director initialization method
  void initializeDirector();
  //GPU initialization method
  void initializeGPU();

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
void LcSim::initializeDirector(){
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

//Initialize memory of the GPU
void LcSim::initializeGPU(){
  //allocate memory for director on GPU
  cudaMalloc( (void**) &gpuDirector
            , 3*iSize_*jSize_*kSize_*sizeof(float));
  //allocate memory for energy on GPU
  cudaMalloc( (void**) &gpuEnergy
            , iSize_*jSize_*kSize_*sizeof(float));
  //allocate memory for torque on GPU
  cudaMalloc( (void**) &gpuTorque
            , 4*3*iSize_*jSize_*kSize_*sizeof(float));
}//initializeGPU
#endif

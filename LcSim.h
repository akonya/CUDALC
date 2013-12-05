#ifndef __LCSIM__
#define __LCSIM__

//  Class for LC simulation
class LcSim{
  //class varriables
  public:
    float *director;
    float *energy;
    int xSize_,ySize_,zSize_;
  //class constructor/destructor
  LcSim(int xSize,int ySize, int zSize);
  ~LcSim();

};

//LcSim constructor
LcSim::LcSim(int xSize, int ySize, int zSize){
  //set sized of cubic grid
  xSize_ = xSize;
  ySize_ = ySize;
  zSize_ = zSize;
  //initialize director and energy arrays
  director = new float[3*xSize_*ySize_*zSize_];
  energy = new float[xSize_*ySize_*zSize_];
}//LcSim

//LcSim deconstructor
LcSim::~LcSim(){
  delete director;
  director = NULL;
  delete energy;
  energy = NULL;  
}//LcSim

#endif

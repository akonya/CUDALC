//===========================================
//               --CUDALC--
//
//  Simulation of LC director field on GPU
//
//  Andrew Konya, 2013
//  Selinger Theory Group
//  Liquid Crystal Institute
//  Kent State University
//===========================================

#include "main.h"

int main(){

  //create LcSim object
  LcSim lcSim = LcSim();

  //intialize director
  lcSim.initDirector();

  //initialize the GPU to use for simulation
  lcSim.initDevice();

  //send host data to the GPU
  lcSim.sendDataToDevice();

  //dynamics loop
  for(int step=0;step<NSTEPS;step++){
    //calculate torque
    lcSim.calculateTorque(0.1);

    //update director
    lcSim.updateDirector();

    //print a frame every FRAMERATE steps
    if(step%FRAMERATE==0){
      //copy data from GPU to CPU
      lcSim.getDataFromDevice();
      //print VTK file 
      lcSim.printVtkFrame(step);
      lcSim.printBmpFrame(step);
    }//if frame
  }//step

  //end LcSim and free memeorye
  lcSim.shutdown();
}//end int main

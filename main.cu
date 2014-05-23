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

  /************************************/
  /* dynamics loop                    */
  /************************************/
  float KbT=0.001;
  int Trate=100000;
  int counter=0;
  for(int step=0;step<NSTEPS;step++){
    //calcualte KbT(step)
//    KbT = 0.01+0.99*exp(-0.00006*float(step));
    //calculate torque
    lcSim.calculateTorque(KbT);

    //update director
    lcSim.updateDirector();

    //print a frame every FRAMERATE steps
    if(step%FRAMERATE==0){

      //copy data from GPU to CPU
      lcSim.getDataFromDevice();

      //print VTK file 
      lcSim.printVtkFrame(step);

      //print BMP file
      lcSim.printBmpFrame(counter);
      
      counter=counter+1;
    }//if frame

    if(step%Trate==0){
      if(KbT<1.0){
        KbT=KbT+0.1;
      }
      if(KbT>=1.0){
        KbT=KbT-0.1;
      }  
    }

  }//step

  //end LcSim and free memeorye
  lcSim.shutdown();
}//end int main

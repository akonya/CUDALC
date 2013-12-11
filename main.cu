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
  LcSim lcSim = LcSim(); //10x10x10 simulation

  //intialize director
  lcSim.initDirector();

  //initialize the GPU to use for simulation
  lcSim.initDevice();

  //send host data to the GPU
  lcSim.sendDataToDevice();

  //do single GPU torque calculation
  lcSim.calculateTorque();

  //do single GPU update
  lcSim.updateDirector(); 

  //get data back from device
  lcSim.getDataFromDevice();
}//end int main

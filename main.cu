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
  LcSim lcSim = LcSim(10,10,10); //10x10x10 simulation

  //initialize the GPU to use for simulation
  lcSim.initializeGPU();

}//end int main

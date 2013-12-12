#ifndef __LCSIMPARAMETERS__
#define __LCSIMPARAMETERS__

#define PI 3.141592654

#define ISIZE 30     //lattice size in i (x) direction
#define JSIZE 30     //lattice size in j (y) direction
#define KSIZE 30     //lattice size in k (z) direction

#define PITCH 10     //pitch in units of lattice spacing
#define QZERO 2.0*PI/PITCH    //2*PI/pitch

#define DELTAT 0.0001    //timestep
#define NSTEPS 10000    //total steps
#define FRAMERATE 1000   //how many steps between frames

#endif

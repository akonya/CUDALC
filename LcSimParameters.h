#ifndef __LCSIMPARAMETERS__
#define __LCSIMPARAMETERS__

#define PI 3.141592654

#define ISIZE 300     //lattice size in i (x) direction
#define JSIZE 300     //lattice size in j (y) direction
#define KSIZE 12     //lattice size in k (z) direction

#define PITCH 12.0     //pitch in units of lattice spacing
#define QZERO 2.0*PI/PITCH    //2*PI/pitch

#define DELTAT 0.0001    //timestep
#define NSTEPS 2000000    //total steps
#define FRAMERATE 10000   //how many steps between frames

//for optics calculations
#define OPTICS_THICKNESS 3.5  //effective optical thickness in microns
#define NeMax           1.7                //max refractive index
#define No              1.5                //ordinary refractive index

#endif

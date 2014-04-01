LC Director Dynamics in CUDA
===============

##Setup on OSC

- clone repository `git clone https://github.com/akonya/CUDALC.git`
- navigate into directory `cd CUDALC`
- load the CUDA module `module load cuda`
- make output folder for VTK files `mkdir VTK`
- make output folder for BMP files `mkdir BMP`
 
You can now run the code with the bach script...

- run using batch script `qsub myscript.job`

OR compile it yourself, and run in a local enviroment

- compile using `nvcc -o runCUDALC main.cu`
- start interactive job `qsub -I -l nodes=1:ppn=1:gpus=1 -l walltime=01:00:00`
- run using `./runCUDALC`

##Primary file contents

- `main.cu` contains the high level logic for the dyanics, ie. main dyanmics loop
- `main.h` contains all the inlcudes for `main.cu`
- `LcSim.h` contains the construction of the `LcSim` class used in the main dynamics loop
- `LcSimParameters.h` contains simulation parameters such as timestep, pitch, size ect.
- `LcSimKernels.h` contains the GPU kernels which are called within the `LcSim` class
- `optics.h`  contains functions to do optics calculations called within `LcSim.printBmpFrame()`




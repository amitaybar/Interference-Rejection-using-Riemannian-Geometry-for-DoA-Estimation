This is the code for the paper "On Interference-Rejection using Riemannian Geometry for Direction of Arrival Estimation", A. Bar and R. Talmon

The code requires the RIR simulator to be installed. Please see more details at https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

The code was tested using Matlab 2018b on Windows 10.

The main file is LoopWrapperMain.m. It calls all the other functions and scripts.
LoopWrapperMain.m calls LoopWrapper.m with an array of flags indicating which figures to plot. The file LoopWrapper.m could be run by itself.

Inside LoopWrapper.m, in lines 25-30, there are flags controlling which figures to plot. Please note that only one of the flags should have a logical '1' value at a time. 



 
  

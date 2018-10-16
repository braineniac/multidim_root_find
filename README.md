A small C program for solving non linear multidimensional systems.  

Prerequisites:   
Install GSL, cmake, clang for your distro.  

Then run:  
mkdir build  
cd build  
cmake ..  
make  

Launch:  
./multidim_root_find [-d] [-o] function  

-o: writed output to a file  
-d: sets if using discreet or standard Newton to find the root  

Available functions:  
five_lin  
five_sq  
five_trig  
five_exp  
rosenbrock  
powell

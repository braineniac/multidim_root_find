A small C program for solving non linear multidimensional systems.

Compilation:
Install GSL, cmake, clang for your distro.

Then run:
mkdir build
cd build
cmake ..
make

./multidim_root_find [-d] [-o] function

-o: writed output to a file
-d: sets if using the Broyden-Newton with the jacobian or just Newton to find the root

Available functions:
five_lin
five_sq
five_trig
five_exp
rosenbrock
powell

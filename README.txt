# hvar
Program for determining frame of minimum Hubble variance
Nov 2014 - Feb 2016
J. H. McKay

The analysis in this program is associated with the work http://arxiv.org/abs/1503.04192 for determining the frame
of minimum Hubble expansion variation.

Quick start:

1) from the base directory run:
$ cd build
$ cmake ..
This may need to more flags to the cmake command if gcc is not the default compiler,
see more detailed instructions in the documentation folder.

2) run
$ make
This will compile the code.

3) Now run either ./main.x or ./Figures.x to see the main results or generate an example of some figures.  If you want to 
choose more figures then open up src/figures.cpp and uncomment the figures you want.  Uncomment all figures to compile all
figures that are possible.  Figure numbers corrspond (roughly) with the numbers in the paper (may need updating +1 to all).



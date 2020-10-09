# Survival-of-the-weakest-in-non-transitive-asymmetric-interactions-among-strains-of-E.-coli
Matlab code for the computational model in the paper


Author : Arianna Miano
Date: June 5, 2020

The code for this manuscript was written in Matlab R2019a (9.6.0.1072779). 
The folder contains Matlab scripts (.m) that can be run on any computer where Matlab R2019a software is installed. The software is available online. 


Once the scripts are opened on Matlab, they can simply be run by clicking the “Run” command in the “Editor” section. 


The scripts used to generate the data for Figure 2 are named "mainRPS_grid_RPS1" and "mainRPS_grid_RPS2", while the remaining files contains the scripts for the parameter sweeps that were used to generate the data in figure 4 and supplementary figure 10.

The code contains different parameters that can be adjusted to simulate different results. For example the initial spatial distribution of “cells” and the relative minimum and maximum probability of death of each cell type. For visualization, if the PLOT variable is set to 1 the code generates a video of the simulated lattice. If the variable SAVE is set to 1, all the lattice frames are saved in the relative folder specified. If the variable VIDEO is set to 1, a video is saved in the same folder at the end of the simulation. 

The density variable can be set to wither 0,1 or 2 to simulate low, medium and high density respectively. 

The variable names for the RPS1 script are as follows:
-P corresponds to Strain B from the RPS-1 ecology in the manuscript
-R corresponds to Strain R from the RPS-1 ecology in the manuscript
-S corresponds to Strain G from the RPS-1 ecology in the manuscript

The variable names for the RPS2 script are as follows:
-P corresponds to Strain G from the RPS-2 ecology in the manuscript
-R corresponds to Strain R from the RPS-2 ecology in the manuscript
-S corresponds to Strain B from the RPS-2 ecology in the manuscript



For the scripts in the main folder the expected run time for a time variable of 3000 is of the order of minutes. 

For the scripts that simulate parameter sweeps the computation time is in the order of hours. It’s adviced to run those simulations overnight. 


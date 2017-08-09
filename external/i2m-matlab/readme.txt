%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Package for analyzing multi-electrode data with a Markov-Ising model.  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This package corresponds to the paper: 

Prediction of spatio-temporal patterns of neural activity from
pairwise correlations.

Olivier Marre, Sami El Boustani, Yves Frégnac and Alain Destexhe

Physical Review Letters, 2009.  

href=http://arxiv.org/abs/0903.0127


The code here allows to reproduce easily the fig 1, analyze your own multi-electrode data, and generate surrogate data with the same statistics than the ones captured by the Markov model. The approach is the following:
-load or generate a raster
-compute the mean activity of each neuron (m), the instantaneous pairwise correlations (C), and the pairwise correlations between time t and time t+1. 
-estimate the h, J and J1 parameters of the model corresponding to the m, C and C1: by an analytical approximation followed by a gradient descent. This might not be enough for a large number of neurons. 
-Estimating the performance of the fit by comparing the prediction, and the empirical estimation, of the ocurrence rate of different spatio-temporal spiking patterns. This is done for different temporal sizes of these patterns. 

The program "BatchOctestGlauber" is performig all these steps. 

1)How to use this program:
-The best is probably to first have a look on the code which reproduces the figure 1. First launch "i2mPath" to set all the sub-directories. Then "BatchOctestGlauber" will do all the analysis (it takes several minutes), and stores the results in the WorkSpace directory. Then launch "Fig1" to draw the figure. 
-To analyze your data, construct a file "spikes.txt" which contain the spike times, with the format explained in /LoadRaster/LoadRaster.m
-the directory InfoTools contains some simple methods to measure the Kullback-Leibler (KL) and the Jensen Shannon (Djs) divergences. 
-/Surrogate/Surrogate.m will generate some surrogate data having the same statistics than the ones captured by the model. 


2) The code is organized in different directories:
-Common: The core of the program. Contains all the functions needed to fit the model to mean and correlations measured from the data. 
-FigurePlot: routines to plot the Figure 1 of the paper
-Figures: directory where automatically generated figures will be stored. 
-Glauber: to simulate the Glauber model
-InfoTools: contains some simple methods to measure the Kullback-Leibler (KL) and the Jensen Shannon (Djs) divergences. 
-LoadRaster: to load and bin a raster
-surrogate: Surrogate.m will generate some surrogate data having the same statistics than the ones captured by the model. 
-Workspace: where the workspace is stored after running one of the batchs. 






 

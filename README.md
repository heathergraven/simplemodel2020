# simplemodel2020
Simple carbon cycle model used in Graven et al. GBC 2020

Coded in Matlab (Tested in versions 2017a and 2020a)

Main file to run the simulations is runSSPsimulations.m

Other files used in the simulation are:
RHS.m - defines system of ODEs
BDprecalc.m - performs initial calculations and loads data
cchems.m, cchems_co2out.m and chemi.m - perform marine carbonate chemistry calculations

PlotOutput.m creates example plots of the simulation output 

Example output file used by PlotOutput.m is SSPsims_200724084931.mat

39 text files have emissions, atmospheric composition and discrimination used as inputs to the model


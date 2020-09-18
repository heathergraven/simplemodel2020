%-----------------------------------------------------------------------
% Return the excess pco2 as a function of cm (DIC in mmol/m**3)
% Heather Graven 2020
%-----------------------------------------------------------------------

function [pco2excess] = chemi(cm,sst,pa0)
     
% borate 
sb=409.07;
% silicate
ssi=46.5;
% phosphate
sp=1.43;
% alkalinity
alk=2333.;
% salinity
sal=35.;
% ah
ah=1.e-8;
     
pco2excess = cchems(cm,sb,ssi,sp,alk,sst,sal,ah) - pa0;



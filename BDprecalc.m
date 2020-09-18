% Perform preliminary calculations to begin carbon cycle simulation
% Heather Graven 2020

global pa0 cm0 cm013 pa013 cm014 pa014 ...
    xi rs R14s akam akma sst0 h1 h2 hm epsab ...
    sal rmbio rmbio14 ... 
    gtc2ppm C14prog prodoc nepsab ntempfb ndifffb ...
    timespan inputdata   cbio01 cbio02 cbio03 ...
    cbio0131 cbio0132 cbio0133 cbio0141 cbio0142 cbio0143 ...
    akba1 akba2 akba3 akab1 akab2 akab3

% set y0 as initial condition with all zeros 
y0=zeros(1,3*49);

% Declare some constants:
gtc2ppm=290/615.6;  % conversion from GtC to ppm
aoc=361419000e6;  % square meter area of the ocean
hm=75.;  % mixed layer depth
h1=25.;  % thermocline layers depth
h2=545.8;  % deep layers depth
rs=0.0112372;  % abundance of C13 vs C12
R14s=1.176e-12;  % abundance of C14 vs C
sal=35;
epsab=-18.;  % land biospheric fractionation (permil)
epsmbio=-18.;  % ocean biospheric fractionation (permil)
xi=10.;  % buffer factor (ignored if nbuff=1)

akma=akam/gtc2ppm/aoc*1.e18/12/hm; % units (1/yr)(mmol/m3)(1/ppm)  
    
% Read in atmospheric data and define initial atmospheric values /
% reference values
[inputdata.CO2time,inputdata.CO2data]=BDreaddata(pco2atmf); timespan.CO2=[inputdata.CO2time(1) inputdata.CO2time(end)]; 
[inputdata.C13time,inputdata.C13data]=BDreaddata(del13atmf); timespan.C13=[inputdata.C13time(1) inputdata.C13time(end)];
[inputdata.C14time,inputdata.C14data]=BDreaddata(del14atmf); timespan.C14=[inputdata.C14time(1) inputdata.C14time(end)];
pa0=inputdata.CO2data(1); 
da0=inputdata.C13data(1); 
Da140=inputdata.C14data(1);
pa013=pa0*rs*(1.+da0/1000.); % 13CO2 approximated using 13C/C rather than 13C/12C
da014 = (Da140 + 2*(da0+25))/(1-2e-3*(da0+25));
pa014 = pa0*(1.+da014/1000.); 

% Read sst data
sst0=18; % set reference value at 18C
if ntempfb==1
    [inputdata.ssttime,inputdata.sstdata]=BDreaddata(sstdatf);
    sst=sst0+inputdata.sstdata(1);
    timespan.sst=[inputdata.ssttime(1) inputdata.ssttime(end)]; 
else
    sst=sst0;
end

% find DIC (mmol/m3) where pCO2 excess is zero: atm in eq with ocean
cm1=1800.;  
cm2=2200.;
tol=.000000001;
optio = optimset('TolX',tol);
[cm0,fval,exitflag]=fzero(@(z) chemi(z,sst,pa0),[cm1 cm2],optio); 

sb=409.07;  % borate
ssi=46.5;  % silicate
sp=1.43;  % phosphate
alk=2333.;  % alkalinity
ah=1.e-8;  % ah
[pm,fco3]=cchems_co3out(cm0,sb,ssi,sp,alk,sst,sal,ah);

% calculate initial/ref fractionation factors and C13 in ocean
% air-sea exchange fractionation factors following Zhang et al 1995
eps_k = -0.86;
eps_aq = +0.0049*sst - 1.31; % in deg C
eps_DIC = 0.014*sst*fco3 - 0.105*sst +10.53;
eps_ao = eps_k + eps_aq; 
eps_oa = eps_k + eps_aq - eps_DIC;
alfaam=eps_ao/1000+1;
alfama=eps_oa/1000+1;
cm013=pa013*(alfaam/alfama)*cm0/pa0;
dm0=(cm013/cm0/rs-1.)*1000.; 

% calculate initial/ref fractionation factors and C14 in ocean
alfaam14=(alfaam-1)*2+1;
alfama14=(alfama-1)*2+1;
Dm014=Da140-50;  % use atm-50 per mil as D14C reference value for surface ocean
dm014=(Dm014 + 2*(dm0+25))/(1-2e-3*(dm0+25));
cm014=cm0*(1.+dm014/1000.);

% calculate fractionation factors and isotopic ratios in marine biota
alfamb =1.+epsmbio/1000.;
alfamb14=1+2*epsmbio/1000; 
rmbio = cm013/cm0*alfamb;
rmbio14 = cm014/cm0*alfamb14;
      
% initial/ref biospheric carbon, C13 and C14 content 
cbio01=pa0*akab1/(akba1*gtc2ppm); % initial biospheric carbon in GtC
cbio02=pa0*akab2/(akba2*gtc2ppm); % initial biospheric carbon in GtC
cbio03=pa0*akab3/(akba3*gtc2ppm); % initial biospheric carbon in GtC
if nepsab==1
    [inputdata.epstime,inputdata.epsdata]=BDreaddata(fepsabf); timespan.eps=[inputdata.epstime(1) inputdata.epstime(end)];
    alfaab=1+inputdata.epsdata(1)/1000;
    alfaab14=(alfaab-1)*2+1;
else
    alfaab=1.+epsab/1000.;
    alfaab14=1.+2*epsab/1000.; 
end
cbio0131=cbio01*alfaab*pa013/pa0; % 13C in biosphere
cbio0132=cbio02*alfaab*pa013/pa0; % 13C in biosphere
cbio0133=cbio03*alfaab*pa013/pa0; % 13C in biosphere
cbio0141=cbio01*alfaab14*pa014/pa0;  % 14C in biosphere
cbio0142=cbio02*alfaab14*pa014/pa0;  % 14C in biosphere
cbio0143=cbio03*alfaab14*pa014/pa0;  % 14C in biosphere

% Read in the rest of the data 
[inputdata.prodtime,inputdata.proddata]=BDreaddata(prodco2f); timespan.prod=[inputdata.prodtime(1) inputdata.prodtime(end)];
[inputdata.prodbiotime,inputdata.prodbiodata]=BDreaddata(prodbiof); timespan.prodbio=[inputdata.prodbiotime(1) inputdata.prodbiotime(end)];
[inputdata.C13fosstime,inputdata.C13fossdata]=BDreaddata(c13foss); timespan.C13foss=[inputdata.C13fosstime(1) inputdata.C13fosstime(end)];
[inputdata.beccstime,inputdata.beccsdata]=BDreaddata(beccsf); timespan.beccs=[inputdata.beccstime(1) inputdata.beccstime(end)]; 
if ndifffb==1
    [inputdata.difftime,inputdata.diffdata]=BDreaddata(diffdatf);
    timespan.diff=[inputdata.difftime(1) inputdata.difftime(end)];
end
if prodoc==1
    [inputdata.prodoctime,inputdata.prodocdata]=BDreaddata(prodocf);
    timespan.prodoc=[inputdata.prodoctime(1) inputdata.prodoctime(end)];
end
if C14prog==1
    [inputdata.ex14time,inputdata.ex14data]=BDreaddata(ex14sourcef);
    timespan.ex14=[inputdata.ex14time(1) inputdata.ex14time(end)]; 
end


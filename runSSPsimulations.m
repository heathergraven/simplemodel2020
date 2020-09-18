global ak akam inputdata C14prog C13prog ndeconv timespan ntempfb epsi ...
    nepsab pa0 pa013 cm014 pa014 gtc2ppm nbuff nfoss ndifffb prodoc ...
    rs R14s  h1 h2 hm cbio0141 cbio0142 cbio0143 fossfac ...
    akba1 akba2 akba3 akab1 akab2 akab3

% This script runs the simple carbon cycle model simulations presented in 
% Graven et al. GBC 2020. A spinup and historical simulation is run first, 
% then 6 different future scenarios are simulated. The model is run for 
% different parameter sets that were selected for their correspondence 
% with observations.
% Heather Graven 2020

% set parameter values
% ak akam epsi akab1 akab2 akab3 akba1 akba2 akba3    
param=[4000	10	0.4	27.1	21.6	210.9	2.4	24.8	580.2;
    4000	10	0.4	33.8	25.6	61.9	2	13.4	212.8;
    4000	11	0.4	21.2	31.5	53.0	2.6	24.1	165.7;
    3000	9	0.4	33.8	25.6	61.9	2	13.4	212.8;
    3000	10	0.4	33.8	25.6	61.9	2	13.4	212.8;
    3000	11	0.4	33.8	25.6	61.9	2	13.4	212.8;
    3000	11	0.4	21.2	31.5	53.0	2.6	24.1	165.7;
    4000	9	0.4	38.1	23.0	107.8	2	15.7	353.2;
    3000	9	0.4	38.1	23.0	107.8	2	15.7	353.2;
    5000	11	0.4	23.2	24.4	105.4	2.5	26.2	299.4;
    5000	10	0	33.8	25.6	61.9	2	13.4	212.8;
    4000	11	0.4	23.2	24.4	105.4	2.5	26.2	299.4;
    4000	10	0	33.8	25.6	61.9	2	13.4	212.8;
    4000	11	0	33.8	25.6	61.9	2	13.4	212.8;
    3000	10	0.4	23.2	24.4	105.4	2.5	26.2	299.4;
    3000	11	0.4	23.2	24.4	105.4	2.5	26.2	299.4;
    3000	9	0	33.8	25.6	61.9	2	13.4	212.8;
    3000	10	0	33.8	25.6	61.9	2	13.4	212.8;
    3000	11	0	33.8	25.6	61.9	2	13.4	212.8;
    4000	9	0	38.1	23.0	107.8	2	15.7	353.2;
    4000	11	0	21.2	31.5	53.0	2.6	24.1	165.7;
    3000	9	0	38.1	23.0	107.8	2	15.7	353.2;
    3000	11	0	21.2	31.5	53.0	2.6	24.1	165.7;
    3000	9	0.4	27.1	21.6	210.9	2.4	24.8	580.2;
    3000	10	0.4	27.1	21.6	210.9	2.4	24.8	580.2;
    5000	11	0	23.2	24.4	105.4	2.5	26.2	299.4;
    4000	11	0	23.2	24.4	105.4	2.5	26.2	299.4;
    4000	10	0	27.1	21.6	210.9	2.4	24.8	580.2;
    3000	10	0	23.2	24.4	105.4	2.5	26.2	299.4;
    3000	11	0	23.2	24.4	105.4	2.5	26.2	299.4;
    3000	9	0	27.1	21.6	210.9	2.4	24.8	580.2;
    3000	10	0	27.1	21.6	210.9	2.4	24.8	580.2;
    3000	11	0.4	21.7	20.1	232.6	2.6	34.7	502.4;
    3000	11	0	21.7	20.1	232.6	2.6	34.7	502.4];

% model setup
ndeconv=1; % ndeconv - 0 prognostic, 1 single deconv
nbuff=1; % nbuff - 0 for constant buffer factor, 1 exact chemistry
nfoss=1; % nfoss - 0 for exponential emissions, 1 for emissions from file
fossfac=1.0; % fossfac - fossil fuel emissions scaling factor
ntempfb=1; % ntempfb - 0 for constant T, else variable T from file
ndifffb=0; % ndifffb - 0 for constant diffusivity, else variable diff from file
nepsab=1; % nepsab - 0 for constant fractionat'n with land biota, else frac from file
prodoc=0; % prodoc - 0 for constant marine biosphere pool, else variable mbio from file

% Define datafiles used in all simulations
del13atmf='c13_cmip6_hist.txt';
del14atmf='c14_cmip6_hist.txt';

% list of future scenario names
SSP={'119','126','245','3B','534','5B'};

% initialize output variable
SSPY=zeros(601,3*49,length(param(:,1)),length(SSP));

% run the model for each parameter set
for t=1:length(param(:,1))
    % specify parameters
    ak=param(t,1); akam=1./param(t,2); epsi=param(t,3);
    akab1=1/param(t,4); akab2=1/param(t,5); akab3=1/param(t,6);
    akba1=1/param(t,7); akba2=1/param(t,8); akba3=1/param(t,9);
    
    for rc=1:6
        % Define datafiles used in different scenarios
        prodco2f=['fossil_SSP' SSP{rc} '.txt'];
        prodbiof=['landuse_SSP' SSP{rc} '.txt'];
        pco2atmf=['co2_SSP' SSP{rc} '.txt'];
        beccsf=['beccs_SSP' SSP{rc} '.txt'];
        fepsabf=['epsab_sj_SSP' SSP{rc} '.txt'];
        sstdatf=['sst_SSP' SSP{rc} '.txt'];
        c13foss='c13foss_constantpost2013.txt';
        
        C14prog=0; % diagnose atm D14C for historical period
        C13prog=0; % diagnose atm d13C for historical period
        
        BDprecalc
        
        % run spinup and historical period to 2005, annual output %%%%%%
        tic;
        [time,Y] = ode23tb(@RHS,-10000:1:2005,y0);
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %calculate steady state total atoms for production term
        atoms14(1)=(Y(11501,1+(3-1)*49)+pa014)/gtc2ppm*1e15/12*R14s*6.02e-3;
        atoms14(2)=(Y(11501,2+(3-1)*49)+cm014)*R14s*hm*aoc/1000*6.02e-3;
        for i=3:39; atoms14(i)=(Y(11501,i+(3-1)*49)+cm014)*R14s*h1*aoc/1000*6.02e-3; end
        for i=40:44; atoms14(i)=(Y(11501,i+(3-1)*49)+cm014)*R14s*h2*aoc/1000*6.02e-3; end
        atoms14(45)=(Y(11501,45+(3-1)*49)+cbio0141)*R14s*1e15/12*6.02e-3;
        atoms14(46)=(Y(11501,46+(3-1)*49)+cbio0142)*R14s*1e15/12*6.02e-3;
        atoms14(47)=(Y(11501,47+(3-1)*49)+cbio0143)*R14s*1e15/12*6.02e-3;
        
        % save number of 14C atoms after spinup for cosmogenic production term
        inputdata.ex14time=1500;
        inputdata.ex14data=sum(atoms14(1:47));
        clear atoms14
        timespan.ex14=[1500 1500];
        
        C14prog=1; % prognostic atm D14C for 2005-2100 period
        C13prog=1; % prognostic atm d13C for 2005-2100 period
        
        % run predictions 2005-2100, annual output %%%%%%%%%%%%%%%%%%%%%
        tic;
        [ftime,fY] = ode23tb(@RHS,2005:1:2100,Y(end,:));
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % combine historical and future simulations in output variable SSPY
        SSPY(1:506,:,t,rc)=Y(11501:end,:);
        SSPY(507:601,:,t,rc)=fY(2:end,:);
        
        clear Y fY
    end
end


% calculate atmospheric values to check simulation
d14C=NaN(size(SSPY)); D14C=d14C; d13C=d14C;
d14C(:,1,:,:)=((SSPY(:,1+(3-1)*49,:,:)+pa014)./(SSPY(:,1,:,:)+pa0)-1)*1000;
d13C(:,1,:,:)=((SSPY(:,1+(2-1)*49,:,:)+pa013)./(SSPY(:,1,:,:)+pa0)/rs-1)*1000;
D14C(:,1,:,:)=d14C(:,1,:,:)-2*(d13C(:,1,:,:)+25).*(1+d14C(:,1,:,:)/1000);

% plot atmospheric CO2
simtime=1500:2100;
figure; plot(simtime,squeeze(SSPY(:,1,:,1))+pa0,'g-')
hold on
plot(simtime,squeeze(SSPY(:,1,:,2))+pa0,'m-')
plot(simtime,squeeze(SSPY(:,1,:,3))+pa0,'b-')
plot(simtime,squeeze(SSPY(:,1,:,4))+pa0,'y-')
plot(simtime,squeeze(SSPY(:,1,:,5))+pa0,'r-')
plot(simtime,squeeze(SSPY(:,1,:,6))+pa0,'k-')

% plot atmospheric d13C
figure; plot(simtime,squeeze(d13C(:,1,:,1)),'g-')
hold on
plot(simtime,squeeze(d13C(:,1,:,2)),'m-')
plot(simtime,squeeze(d13C(:,1,:,3)),'b-')
plot(simtime,squeeze(d13C(:,1,:,4)),'y-')
plot(simtime,squeeze(d13C(:,1,:,5)),'r-')
plot(simtime,squeeze(d13C(:,1,:,6)),'k-')
plot(inputdata.C13time,inputdata.C13data,'c-')

% plot atmospheric D14C
figure; plot(simtime,squeeze(D14C(:,1,:,1)),'g-')
hold on
plot(simtime,squeeze(D14C(:,1,:,2)),'m-')
plot(simtime,squeeze(D14C(:,1,:,3)),'b-')
plot(simtime,squeeze(D14C(:,1,:,4)),'y-')
plot(simtime,squeeze(D14C(:,1,:,5)),'r-')
plot(simtime,squeeze(D14C(:,1,:,6)),'k-')
plot(inputdata.C14time,inputdata.C14data,'c-')

% save output in matlab workspace
save(['SSPsims_' datestr(now,'yymmddHHMMSS')])



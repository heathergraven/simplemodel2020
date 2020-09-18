function a = RHS(time, y)

global pa0 cm0  cm013 pa013 cm014 pa014 ...
    xi rs R14s  epsi akam akma sst0 h1 h2 hm ak epsab ...
    sal rmbio rmbio14 ... 
    nfoss fossfac gtc2ppm nbuff C14prog C13prog nepsab ndeconv ntempfb ndifffb prodoc ...
    timespan inputdata cbio01 cbio02 cbio03 ...
    cbio0131 cbio0132 cbio0133 cbio0141 cbio0142 cbio0143 ...
    akba1 akba2 akba3 akab1 akab2 akab3

% This script defines the ordinary differential equations in the model
% Heather Graven 2020

% initialize array of ode's
a=zeros(49*3,1); 
% 1: atmosphere
% 2: mixed layer
% 3-39: thermocline
% 40-44: deep ocean
% 45-47: biosphere boxes
% 48-49: unused
% First 49 correspond to carbon, next 49 to 13C and last 49 to 14C

% set temperature
if ntempfb==1  % if variable sst
    if time<timespan.sst(1)
        sstt=inputdata.sstdata(1);
    elseif time>=timespan.sst(1) && time<=timespan.sst(2)
        sstt=interp1(inputdata.ssttime,inputdata.sstdata,time);
    else
        sstt=inputdata.sstdata(end);
    end
    ssts=sstt+sst0;
else  % if constant sst
    ssts=sst0;
end

% set carbonate chemistry
if nbuff==0  % if constant revelle factor
    pm=pa0*(1.+xi*y(2)/cm0);  %y(2,1) 
else  % if chemistry calculated explicitly     
    sb=409.07;  % borate
    ssi=46.5;  % silicate
    sp=1.43;  % phosphate
    alk=2333.;  % alkalinity
    ah=1.e-8;  % ah
    [pm,fco3]=cchems_co3out(cm0+y(2),sb,ssi,sp,alk,ssts,sal,ah);
end

% calculate air-sea exchange fractionation factors following Zhang et al 1995
eps_k = -0.86;
eps_aq = +0.0049*ssts - 1.31; % in deg C
eps_DIC = 0.014*ssts*fco3 - 0.105*ssts +10.53;
eps_ao = eps_k + eps_aq; 
eps_oa = eps_k + eps_aq - eps_DIC;
alfaam=eps_ao/1000+1;
alfama=eps_oa/1000+1;
alfaam14=(alfaam-1)*2+1;
alfama14=(alfama-1)*2+1;
  
% calculate biospheric ratios to use in land use terms
rbio1=(cbio0131+y(45+(2-1)*49))/(cbio01+y(45));  
rbio2=(cbio0132+y(46+(2-1)*49))/(cbio02+y(46));  
rbio3=(cbio0133+y(47+(2-1)*49))/(cbio03+y(47));  
rbio141=(cbio0141+y(45+(3-1)*49))/(cbio01+y(45));  
rbio142=(cbio0142+y(46+(3-1)*49))/(cbio02+y(46));  
rbio143=(cbio0143+y(47+(3-1)*49))/(cbio03+y(47));  

% set fossil fuel emissions
if nfoss==0  % if exponential emissions
    amu=1./34.; p0=5.2;    
    q=fossfac*p0*exp(amu*(time-1980.))*gtc2ppm;
elseif nfoss==1  % if reported emissions
    if time<timespan.prod(1)
        q=fossfac*inputdata.proddata(1)*gtc2ppm;
    elseif time>=timespan.prod(1) && time<=timespan.prod(2)
        q=fossfac*interp1(inputdata.prodtime,inputdata.proddata,time)*gtc2ppm;
    else
        q=fossfac*inputdata.proddata(end)*gtc2ppm;
    end
end

% set beccs
if time<timespan.beccs(1)
    beccs=inputdata.beccsdata(1)*gtc2ppm;
elseif time>=timespan.beccs(1) && time<=timespan.beccs(2)
    beccs=interp1(inputdata.beccstime,inputdata.beccsdata,time)*gtc2ppm;
else
    beccs=inputdata.beccsdata(end)*gtc2ppm;
end

% set land use emissions
if time<timespan.prodbio(1)
    qbio=inputdata.prodbiodata(1)*gtc2ppm;
elseif time>=timespan.prodbio(1) && time<=timespan.prodbio(2)
    qbio=interp1(inputdata.prodbiotime,inputdata.prodbiodata,time)*gtc2ppm;
else
    qbio=inputdata.prodbiodata(end)*gtc2ppm;
end
% attribute fraction of reported ff emissions to land use emissions
if ndeconv==0
    qbio=qbio+(1-fossfac)*q; 
end  

% set marine biosphere change
if prodoc==1
    qmbio=interp1(inputdata.prodoctime,inputdata.prodocdata,time);
else
    qmbio=0; 
end

% set air-land exchange fractionation factors
if nepsab==1  % if variable fractionation
    if time<timespan.eps(1)
        alfaab =1.+inputdata.epsdata(1)/1000.;
        alfaab14=1.+2*inputdata.epsdata(1)/1000.;
    elseif time>=timespan.eps(1) && time<=timespan.eps(2)
        alfaab =1.+interp1(inputdata.epstime,inputdata.epsdata,time)/1000.;
        alfaab14=1.+2*interp1(inputdata.epstime,inputdata.epsdata,time)/1000.;
    else
        alfaab =1.+inputdata.epsdata(end)/1000.;
        alfaab14=1.+2*inputdata.epsdata(end)/1000.;
    end
else  % if constant fractionation
    alfaab=1.+epsab/1000.;
    alfaab14=1.+2*epsab/1000.; 
end

% calculate derivatives by linear interpolation with 1 year window
dt=0.5;
% read observed atmospheric CO2 and derivatives
if time<timespan.CO2(1)+dt
    xpa=0;
    dxpadt=0;
elseif time>=timespan.CO2(1)+dt && time<=timespan.CO2(2)-dt
    xpa=interp1(inputdata.CO2time,inputdata.CO2data,time)-pa0;
    xpa1=interp1(inputdata.CO2time,inputdata.CO2data,time-dt)-pa0;
    xpa2=interp1(inputdata.CO2time,inputdata.CO2data,time+dt)-pa0;
    dxpadt=(xpa2-xpa1)/(2*dt);
else
    xpa=inputdata.CO2data(end)-pa0;
    dxpadt=0;
end

% read observed atmospheric d13C and derivative
if time<timespan.CO2(1)+dt
    xda=0;
    ddadt=0;
elseif time>=timespan.CO2(1)+dt && time<timespan.C13(1)+dt
    xda=pa013/pa0*(xpa+pa0)-pa013;
    xda1=pa013/pa0*(xpa1+pa0)-pa013;
    xda2=pa013/pa0*(xpa2+pa0)-pa013;
    ddadt=(xda2-xda1)/(2*dt);
elseif time>=timespan.C13(1)+dt && time<=timespan.CO2(2)-dt
    xda=(xpa+pa0)*rs*(1+interp1(inputdata.C13time,inputdata.C13data,time)/1000)-pa013;
    xda1=(xpa1+pa0)*rs*(1+interp1(inputdata.C13time,inputdata.C13data,time-dt)/1000)-pa013;
    xda2=(xpa2+pa0)*rs*(1+interp1(inputdata.C13time,inputdata.C13data,time+dt)/1000)-pa013;
    ddadt=(xda2-xda1)/(2*dt);
else
    xda=(xpa+pa0)*rs*(1+inputdata.C13data(end)/1000)-pa013;
    ddadt=0;
end

% read observed atmospheric D14C and derivative
if time<timespan.C14(1)+dt
    xda14=0;
    ddadt14=0;
elseif time>=timespan.C14(1)+dt && time<timespan.CO2(1)+dt
    xda14=(((interp1(inputdata.C14time,inputdata.C14data,time)+...
        2*((pa013/pa0/rs-1)*1000+25))/(1-2e-3*...
        ((pa013/pa0/rs-1)*1000+25)))/1000+1)*pa0...
        -pa014;
    ddadt14=((((interp1(inputdata.C14time,inputdata.C14data,time+dt)+...
        2*((pa013/pa0/rs-1)*1000+25))/(1-2e-3*...
        ((pa013/pa0/rs-1)*1000+25)))/1000+1)*pa0-...
        (((interp1(inputdata.C14time,inputdata.C14data,time-dt)+...
        2*((pa013/pa0/rs-1)*1000+25))/(1-2e-3*...
        ((pa013/pa0/rs-1)*1000+25)))/1000+1)*pa0)/(2*dt);
elseif time>=timespan.CO2(1)+dt && time<=timespan.C14(2)-dt
    xda14=(((interp1(inputdata.C14time,inputdata.C14data,time)+...
        2*(((xda+pa013)/(xpa+pa0)/rs-1)*1000+25))/(1-2e-3*...
        (((xda+pa013)/(xpa+pa0)/rs-1)*1000+25)))/1000+1)*(xpa+pa0)...
        -pa014;
    xda141=(((interp1(inputdata.C14time,inputdata.C14data,time-dt)+...
        2*(((xda1+pa013)/(xpa1+pa0)/rs-1)*1000+25))/(1-2e-3*...
        (((xda1+pa013)/(xpa1+pa0)/rs-1)*1000+25)))/1000+1)*(xpa1+pa0)...
        -pa014;
    xda142=(((interp1(inputdata.C14time,inputdata.C14data,time+dt)+...
        2*(((xda2+pa013)/(xpa2+pa0)/rs-1)*1000+25))/(1-2e-3*...
        (((xda2+pa013)/(xpa2+pa0)/rs-1)*1000+25)))/1000+1)*(xpa2+pa0)...
        -pa014;
    ddadt14=(xda142-xda141)/(2*dt);
else
    xda14=(inputdata.C14data(end)+...
        2*(((xda+pa013)/(xpa+pa0)/rs-1)*1000+25))/(1-2e-3*...
        (((xda+pa013)/(xpa+pa0)/rs-1)*1000+25))-pa014;
    ddadt14=0;
end 

% read 14C atom data
if C14prog==1
    if time<timespan.ex14(1)+dt
        ex14c=inputdata.ex14data(1);
        dex14cdt=0;
    elseif time>=timespan.ex14(1)+dt && time<=timespan.ex14(2)-dt
        ex14c=interp1(inputdata.ex14time,inputdata.ex14data,time);
        dex14cdt=(interp1(inputdata.ex14time,inputdata.ex14data,time+dt)-...
        interp1(inputdata.ex14time,inputdata.ex14data,time-dt))/(2*dt);
    else
        ex14c=inputdata.ex14data(end);
        dex14cdt=0;
    end
end
   
% set 13C/12C ratio in fossil fuel carbon
if nfoss==0        
    qf=-25;
    rfosstime=rs*(1+qf/1000.);
elseif nfoss==1    
    if time<timespan.C13foss(1)
        rfosstime=rs*(1+inputdata.C13fossdata(end)/1000.);
    elseif time>=timespan.C13foss(1) && time<=timespan.C13foss(2)
        rfosstime=rs*(1+interp1(inputdata.C13fosstime,inputdata.C13fossdata,time)/1000.);
    else
        rfosstime=rs*(1+inputdata.C13fossdata(end)/1000.);
    end
end

% constants of discretization scheme in ocean
if ndifffb==1
    if time<timespan.diff(1)
        aks=inputdata.diffdata(1);
    elseif time>=timespan.diff(1) && time<=timespan.diff(2)
        aks=interp1(inputdata.difftime,inputdata.diffdata,time);
    else
        aks=inputdata.diffdata(end);
    end
else
    aks=ak;
end
ak1=aks/(h1*h1);
ak2=aks/(h2*h2);
akv=2*aks/(h1*(h1+h2));
akn=2*aks/(h2*(h1+h2));
akmd=2*aks/(hm*(hm+h1));
akdm=2*aks/(h1*(hm+h1));


% now program the equations
if ndeconv==0   % Forward Model run
    % a(1,1) = a(1)  Atmospheric CO2
    a(1)=akam*(pm-(pa0+y(1)))+q+beccs+qbio-...
        akab1*(pa0+epsi*y(1))+akba1*(y(45)+cbio01)*gtc2ppm- ...
        akab2*(pa0+epsi*y(1))+akba2*(y(46)+cbio02)*gtc2ppm- ...
        akab3*(pa0+epsi*y(1))+akba3*(y(47)+cbio03)*gtc2ppm;  

    % a(1,2) = a(50)  Atmospheric 13C
    a(1+(2-1)*49)=akam*(alfama*((y(2+(2-1)*49)+cm013)/(y(2)+cm0))*pm-...
        alfaam*(y(1+(2-1)*49)+pa013))+...
        (akba1*(y(45+(2-1)*49)+cbio0131)*gtc2ppm-alfaab*akab1*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        (akba2*(y(46+(2-1)*49)+cbio0132)*gtc2ppm-alfaab*akab2*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        (akba3*(y(47+(2-1)*49)+cbio0133)*gtc2ppm-alfaab*akab3*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        (q*rfosstime) ...
        +beccs*(alfaab*(pa013+y(1+(2-1)*49))/(pa0+y(1))) ...
        +(qbio>=0)*(qbio*rbio2) ...
        +(qbio<0)*(alfaab*qbio*(pa013+y(1+(2-1)*49))/(pa0+y(1)));
    
    % Atmospheric 14C
    if C14prog==0
        a(1+(3-1)*49)=ddadt14;
    elseif C14prog==1
        a(1+(3-1)*49)=akam*(alfama14*((y(2+(3-1)*49)+cm014)/(y(2)+cm0))*pm...
            -alfaam14*(y(1+(3-1)*49)+pa014) )...
            +(akba1*(y(45+(3-1)*49)+cbio0141)*gtc2ppm-...
            alfaab14*akab1*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            +(akba2*(y(46+(3-1)*49)+cbio0142)*gtc2ppm-...
            alfaab14*akab2*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            +(akba3*(y(47+(3-1)*49)+cbio0143)*gtc2ppm-...
            alfaab14*akab3*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            -(y(1+(3-1)*49)+pa014)/8267....
            +(qbio>=0)*(qbio*rbio142)...
            +(qbio<0)*(alfaab14*qbio*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
            +beccs*(alfaab14*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
            +ex14c*gtc2ppm/1e15*12/R14s/6.02e-3/8267 ...
            +0.30852866; % assume nuclear production fixed at 2008 value
        % assume negative fossil fuel emissions result from bio-CCS
    end
    
    % Biosphere
    % a(45,1) = a(45)  Land Biosphere Total Carbon - Box 1
    a(45)=akab1*(pa0+epsi*y(1))/gtc2ppm-akba1*(y(45)+cbio01);
    % a(46,1) = a(46)  Land Biosphere Total Carbon - Box 2 - Land use
    % fluxes go here
    a(46)=akab2*(pa0+epsi*y(1))/gtc2ppm-akba2*(y(46)+cbio02)-qbio/gtc2ppm;
    % a(47,1) = a(47)  Land Biosphere Total Carbon - Box 3
    a(47)=akab3*(pa0+epsi*y(1))/gtc2ppm-akba3*(y(47)+cbio03);
    
    
    % a(45,2) = a(94)  Land Biosphere 13C - box 1
    a(45+(2-1)*49)= -(akba1*(y(45+(2-1)*49)+cbio0131)-...
        alfaab*akab1*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm);
    % a(46,2)  Land Biosphere 13C - box 2 - land use fluxes go
    % here
    a(46+(2-1)*49)= -(akba2*(y(46+(2-1)*49)+cbio0132)-...
        alfaab*akab2*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(qbio>=0)*(qbio/gtc2ppm*rbio2) ...
        -(qbio<0)*(alfaab*qbio/gtc2ppm*(pa013+y(1+(2-1)*49))/(pa0+y(1)));
    % a(47,2)  Land Biosphere 13C - box 3
    a(47+(2-1)*49)= -(akba3*(y(47+(2-1)*49)+cbio0133)-...
        alfaab*akab3*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm);
  
    %  Land Biosphere 14C - box 1
    a(45+(3-1)*49)= -(akba1*(y(45+(3-1)*49)+cbio0141)-...
        alfaab14*akab1*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(y(45+(3-1)*49)+cbio0141)/8267.;
    %  Land Biosphere 14C - box 2
    a(46+(3-1)*49)= -(akba2*(y(46+(3-1)*49)+cbio0142)-...
        alfaab14*akab2*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(qbio>=0)*(qbio/gtc2ppm*rbio142)...
        -(qbio<0)*(alfaab14*qbio*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
        -(y(46+(3-1)*49)+cbio0142)/8267.;
    %  Land Biosphere 14C - box 3
    a(47+(3-1)*49)= -(akba3*(y(47+(3-1)*49)+cbio0143)-...
        alfaab14*akab3*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(y(47+(3-1)*49)+cbio0143)/8267.;


    % Mixed Layer
    % a(2,1) = a(2)  Mixed Layer total carbon
    a(2)=-akma*(pm-(pa0+y(1)))+akmd*(-3*y(2)+4*y(3)-y(4))+akma/akam*qmbio*gtc2ppm;

    % a(2,2) = a(51)  Mixed Layer 13C
    a(2+(2-1)*49)=-akma*(alfama*((y(2+(2-1)*49)+cm013)/(y(2)+cm0))*pm-...
        alfaam*(y(1+(2-1)*49)+pa013))+akmd*(-3*y(2+(2-1)*49)+4*y(3+(2-1)*49)-...
        y(4+(2-1)*49))+akma/akam*qmbio*gtc2ppm*rmbio;

    % Mixed Layer 14C
    a(2+(3-1)*49)=-akma*(alfama14*((y(2+(3-1)*49)+cm014)/(y(2)+cm0))*pm-...
        alfaam14*(y(1+(3-1)*49)+pa014))+akmd*(-3*y(2+(3-1)*49)+4*y(3+(3-1)*49)-...
        y(4+(3-1)*49))+akma/akam*qmbio*gtc2ppm*rmbio14...
        -(y(2+(3-1)*49)+cm014)/8267.;
      
elseif ndeconv==1  % Single Deconvolution
    % a(1,1) = a(1)  Atmospheric CO2
    a(1)=dxpadt; 
    sdbioflux=(akam*(pm-(pa0+y(1)))+q+beccs+qbio- ...
        akab1*(pa0+epsi*y(1))+akba1*(y(45)+cbio01)*gtc2ppm- ...
        akab2*(pa0+epsi*y(1))+akba2*(y(46)+cbio02)*gtc2ppm- ...
        akab3*(pa0+epsi*y(1))+akba3*(y(47)+cbio03)*gtc2ppm)-dxpadt;
    % assuming extra flux goes into or out of biosphere
    
    % a(45,1) = a(45)  Land Biosphere Total Carbon - Box 1
    a(45)=akab1*(pa0+epsi*y(1))/gtc2ppm-akba1*(y(45)+cbio01);
    
    % a(46,1) = a(46)  Land Biosphere Total Carbon - Box 2 - Land use and
    % residual fluxes go here
    a(46)=akab2*(pa0+epsi*y(1))/gtc2ppm-akba2*(y(46)+cbio02)-qbio/gtc2ppm+sdbioflux/gtc2ppm;
    
    % a(47,1) = a(47)  Land Biosphere Total Carbon - Box 3
    a(47)=akab3*(pa0+epsi*y(1))/gtc2ppm-akba3*(y(47)+cbio03);
    
    if sdbioflux>0 % flux is into biosphere
        sdbio13flux=sdbioflux*alfaab*(pa013+y(1+(2-1)*49))/(pa0+y(1));
        sdbio14flux=sdbioflux*alfaab14*(pa014+y(1+(3-1)*49))/(pa0+y(1));
    else % flux is out of biosphere
        sdbio13flux=sdbioflux*rbio2;
        sdbio14flux=sdbioflux*rbio142;
    end
    
    if C13prog==0
        a(1+(2-1)*49)=ddadt;
    else
    % atmospheric d13C 
    a(1+(2-1)*49)=akam*(alfama*((y(2+(2-1)*49)+cm013)/(y(2)+cm0))*pm-...
        alfaam*(y(1+(2-1)*49)+pa013))+...
        (akba1*(y(45+(2-1)*49)+cbio0131)*gtc2ppm-alfaab*akab1*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        (akba2*(y(46+(2-1)*49)+cbio0132)*gtc2ppm-alfaab*akab2*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        (akba3*(y(47+(2-1)*49)+cbio0133)*gtc2ppm-alfaab*akab3*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1)))+...
        -sdbio13flux ...
        +(q*rfosstime) ...
        +beccs*(alfaab*(pa013+y(1+(2-1)*49))/(pa0+y(1))) ...
        +(qbio>=0)*(qbio*rbio2) ...
        +(qbio<0)*(alfaab*qbio*(pa013+y(1+(2-1)*49))/(pa0+y(1)));
    end

    % Atmospheric 14C
    if C14prog==0
        a(1+(3-1)*49)=ddadt14;
    elseif C14prog==1
        a(1+(3-1)*49)=akam*(alfama14*((y(2+(3-1)*49)+cm014)/(y(2)+cm0))*pm...
            -alfaam14*(y(1+(3-1)*49)+pa014) )...
            +(akba1*(y(45+(3-1)*49)+cbio0141)*gtc2ppm-...
            alfaab14*akab1*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            +(akba2*(y(46+(3-1)*49)+cbio0142)*gtc2ppm-...
            alfaab14*akab2*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            +(akba3*(y(47+(3-1)*49)+cbio0143)*gtc2ppm-...
            alfaab14*akab3*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1)))...
            -(y(1+(3-1)*49)+pa014)/8267....
            -sdbio14flux ...
            +(qbio>=0)*(qbio*rbio142)...
            +(qbio<0)*(alfaab14*qbio*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
            +beccs*(alfaab14*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
            +ex14c*gtc2ppm/1e15*12/R14s/6.02e-3/8267 ...
            +0.30852866; % assume nuclear production fixed at 2008 value
        % assume negative fossil fuel emissions result from bio-CCS
    end

    % a(45,2) = a(94)  Land Biosphere 13C - box 1
    a(45+(2-1)*49)= -(akba1*(y(45+(2-1)*49)+cbio0131)-...
        alfaab*akab1*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm);
    % a(46,2)  Land Biosphere 13C - box 2 - land use and residual flux go
    % here
    a(46+(2-1)*49)= -(akba2*(y(46+(2-1)*49)+cbio0132)-...
        alfaab*akab2*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        +sdbio13flux/gtc2ppm ...
        -(qbio>=0)*(qbio/gtc2ppm*rbio2) ...
        -(qbio<0)*(alfaab*qbio/gtc2ppm*(pa013+y(1+(2-1)*49))/(pa0+y(1)));
    % a(47,2)  Land Biosphere 13C - box 3
    a(47+(2-1)*49)= -(akba3*(y(47+(2-1)*49)+cbio0133)-...
        alfaab*akab3*(pa013+y(1+(2-1)*49))*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm);


    %  Land Biosphere 14C - box 1
    a(45+(3-1)*49)= -(akba1*(y(45+(3-1)*49)+cbio0141)-...
        alfaab14*akab1*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(y(45+(3-1)*49)+cbio0141)/8267.;
    %  Land Biosphere 14C - box 2
    a(46+(3-1)*49)= -(akba2*(y(46+(3-1)*49)+cbio0142)-...
        alfaab14*akab2*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        +sdbio14flux/gtc2ppm ...
        -(qbio>=0)*(qbio/gtc2ppm*rbio142)...
        -(qbio<0)*(alfaab14*qbio*(y(1+(3-1)*49)+pa014)/(pa0+y(1)))...
        -(y(46+(3-1)*49)+cbio0142)/8267.;
    %  Land Biosphere 14C - box 3
    a(47+(3-1)*49)= -(akba3*(y(47+(3-1)*49)+cbio0143)-...
        alfaab14*akab3*(y(1+(3-1)*49)+pa014)*(pa0+epsi*y(1))/(pa0+y(1))/gtc2ppm)...
        -(y(47+(3-1)*49)+cbio0143)/8267.;

    
    % a(2,1) = a(2)  Mixed Layer total carbon
    a(2)=-akma*(pm-(pa0+y(1)))+akmd*(-3*y(2)+4*y(3)-y(4))+akma/akam*qmbio*gtc2ppm;

    % a(2,2) = a(51)  Mixed Layer 13C
    a(2+(2-1)*49)=-akma*(alfama*((y(2+(2-1)*49)+cm013)/(y(2)+cm0))*pm-...
        alfaam*(y(1+(2-1)*49)+pa013))+akmd*(-3*y(2+(2-1)*49)+4*y(3+(2-1)*49)-...
        y(4+(2-1)*49))+akma/akam*qmbio*gtc2ppm*rmbio;

    % Mixed Layer 14C
    a(2+(3-1)*49)=-akma*(alfama14*((y(2+(3-1)*49)+cm014)/(y(2)+cm0))*pm-...
        alfaam14*(y(1+(3-1)*49)+pa014))+akmd*(-3*y(2+(3-1)*49)+4*y(3+(3-1)*49)-...
        y(4+(3-1)*49))+akma/akam*qmbio*gtc2ppm*rmbio14...
        -(y(2+(3-1)*49)+cm014)/8267.;   
end 

% first thermocline box
a(3)=-akdm*(-3*y(2)+4*y(3)-y(4))+ak1*(y(4)-y(3));
a(3+(2-1)*49)=-akdm*(-3*y(2+(2-1)*49)+4*y(3+(2-1)*49)-...
    y(4+(2-1)*49))+ak1*(y(4+(2-1)*49)-y(3+(2-1)*49));
a(3+(3-1)*49)=-akdm*(-3*y(2+(3-1)*49)+4*y(3+(3-1)*49)-...
    y(4+(3-1)*49))+ak1*(y(4+(3-1)*49)-y(3+(3-1)*49))-(y(3+(3-1)*49)+cm014)/8267.;

% Thermocline   
for i=4:38 
    a(i)=ak1*(y(i-1)-2*y(i)+y(i+1));  % Total Carbon
    a(i+(2-1)*49)=ak1*(y(i-1+(2-1)*49)-2*y(i+(2-1)*49)+y(i+1+(2-1)*49));  % 13C
	a(i+(3-1)*49)=ak1*(y(i-1+(3-1)*49)-2*y(i+(3-1)*49)+y(i+1+(3-1)*49))-(y(i+(3-1)*49)+cm014)/8267.;  % 14C
end
a(39)=ak1*(y(38)-y(39))+akv*(y(40)-y(39));
a(39+(2-1)*49)=ak1*(y(38+(2-1)*49)-y(39+(2-1)*49))+akv*(y(40+(2-1)*49)-y(39+(2-1)*49));
a(39+(3-1)*49)=ak1*(y(38+(3-1)*49)-y(39+(3-1)*49))+akv*(y(40+(3-1)*49)-y(39+(3-1)*49))-(y(39+(3-1)*49)+cm014)/8267.;

% Deep Sea
a(40)=akn*(y(39)-y(40))+ak2*(y(41)-y(40));  % Total Carbon
a(40+(2-1)*49)=akn*(y(39+(2-1)*49)-y(40+(2-1)*49))+ak2*(y(41+(2-1)*49)-y(40+(2-1)*49));  % 13C
a(40+(3-1)*49)=akn*(y(39+(3-1)*49)-y(40+(3-1)*49))+ak2*(y(41+(3-1)*49)-y(40+(3-1)*49))-(y(40+(3-1)*49)+cm014)/8267.;  % 14C
      
for i=41:43
    a(i)=ak2*(y(i-1)-2*y(i)+y(i+1));  % Total Carbon
    a(i+(2-1)*49)=ak2*(y(i-1+(2-1)*49)-2*y(i+(2-1)*49)+y(i+1+(2-1)*49));  % 13C
	a(i+(3-1)*49)=ak2*(y(i-1+(3-1)*49)-2*y(i+(3-1)*49)+y(i+1+(3-1)*49))-(y(i+(3-1)*49)+cm014)/8267.;  % 14C
end 

a(44)=ak2*(y(43)-y(44));  % Total Carbon
a(44+(2-1)*49)=ak2*(y(43+(2-1)*49)-y(44+(2-1)*49));  % 13C
a(44+(3-1)*49)=ak2*(y(43+(3-1)*49)-y(44+(3-1)*49))-(y(44+(3-1)*49)+cm014)/8267.;  % 14C
     



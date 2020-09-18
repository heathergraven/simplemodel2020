%------------------------------------------------------
% Carbonate chemistry calculation 
%
% input:
%
% sc	total carbon
% sb	total borate
% ssi	total silicium
% sp	total phosphorus
% alk	total alkalinity
% t     temperature (degree Celsius)
% sal   salinity (permil)
% ah	first guess activity of hydrogen (must be greater than zero)
%
% output:
%
% ah	resulting hydrogen activity
% ph	pH
% co2	
% hco3
% co3
% pco2
%
% all concentrations are given in umol/kg (mmol/m**3).
%
% chemical constants as stated in Peng et al., 1987
%
% Heather Graven 2020
%------------------------------------------------------

function [pco2,fco3] = cchems_co3out(sc,sb,ssi,sp,alk,t,sal,ah)

nmax=1000000000000;
eps=.00000001;

ta=t+273.15;
tac=ta/100;

% Reaction constants calculated for given Temp and Salinity
ak1=10^(13.7201-0.031334*ta-3235.76/ta-1.3e-5*sal*ta+0.1032*sqrt(sal));
a1=-5371.9645-1.671221*ta+128375.28/ta;
a2=2194.3055*log10(ta)-0.22913*sal-18.3802*log10(sal);
a3=8.0944e-4*sal*ta+5617.11*log10(sal)/ta-2.136*sal/ta;
ak2=10^(a1+a2+a3);

akb=10^(-9.26+0.00886*sal+0.01*t);

aksi=4*10^-10;

akp1=2*10^-2;
akp2=exp(-9.039-1450/ta);
akp3=exp(4.466-7276/ta);

akw=exp(148.9802-13847.26/ta-23.6521*log(ta)-0.019813*sal+sqrt(sal)*(-79.2447+3298.72/ta+12.0408*log(ta)));
fh=1.29-0.00204*ta+sal*sal*(4.61e-4-1.48e-6*ta);
      
alfa=exp(-60.2409+9345.17/ta+23.3585*log(tac)+sal*(0.023517-0.023656*tac+0.0047036*tac*tac));
      
n=0;
a=1/ah;
sum(1)=sc;
sum(2)=sb;
sum(3)=ssi;
sum(4)=sp;
alka=alk;

x=0;
while x==0
    n=n+1;
    al=a;
   
    g(1,1)=1+ak1*a*(1+ak2*a);
    g(1,2)=ak1*(1+2*ak2*a);
    g(1,3)=2*ak1*ak2;

    g(2,1)=1+akb*a;
    g(2,2)=akb;
    g(2,3)=0;

    g(3,1)=1+aksi*a;
    g(3,2)=aksi;
    g(3,3)=0;

    g(4,1)=1+akp1*a*(1+akp2*a*(1+akp3*a));
    g(4,2)=akp1*(1+akp2*a*(2+3*akp3*a));
    g(4,3)=akp1*akp2*(2+6*akp3*a);

    h=0;
    hs=0;

    for w=1:4
        h=h+sum(w)*g(w,2)/g(w,1);
        hs=hs+sum(w)*(g(w,3)*g(w,1)-g(w,2)*g(w,2))/(g(w,1)*g(w,1));
    end        

    hs=h+a*hs+1*10^6*akw*fh+1*10^6/(a*a*fh);
    h=a*(h+1*10^6*akw*fh)-1*10^6/(a*fh)-alka;
    a=a-h/hs;
    
    if (n>nmax)
        x=1;
        disp(sprintf('No convergence in module CCHEMS, ah0= %d', ah))
    elseif (abs((a-al)/al)<eps)
        x=1;
    
        ah=1/a;
        hco3=sc/(1.+ak2/ah+ah/ak1);
        co3=hco3*ak2/ah;
        co2=hco3*ah/ak1;
        pco2=co2/alfa;
        ph=-log10(ah/fh);
        
        fco3=co3/(co3+hco3+co2);
    end
    % otherwise, x is still zero and loop restarts
end




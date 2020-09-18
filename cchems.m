%------------------------------------------------------
% Full chemistry calculation using single precision
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

function [pco2] = cchems(sc,sb,ssi,sp,alk,t,sal,ah)

%disp(sprintf('Function CCHEMS started'))

nmax=1000000000000;
eps=.00000001;

ta=t+273.15;
tac=ta/100;


% HG - Reaction constants calculated for given Temp and Salinity
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

   
% this is where the loop starts

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

        %disp(sprintf('Calculated carbonate concentrations \n ph = %d \n pCO2 = %d \n CO2 = %d \n HCO3 = %d \n CO3 = %d', ph, pco2, co2, hco3, co3))
    end

    % HG - otherwise, x is still zero and loop restarts

end




% Original Fortran Code,  Lines 313-445,  bdif7.f
%{
      subroutine cchems(sc,sb,ssi,sp,alk,t,s,ah,
     .                  ph,co2,hco3,co3,pco2)
c-------------------------------------------------------------------     
c full chemistry calculation using single precision
c
c input:
c
c sc	total carbon
c sb	total borate
c ssi	total silicium
c sp	total phosphorus
c alk	total alkalinity
c t	temperature (degree Celsius)
c s	salinity (permil)
c ah	first guess activity of hydrogen (must be greater than zero)
c
c output:
c
c ah	resulting hydrogen activity
c ph	pH
c co2	
c hco3
c co3
c pco2
c
c all concentrations are given in umol/kg (mmol/m**3).
c
c chemical constants as stated in Peng et al., 1987
c
c single precision
c 
c						mh, 11-jan-1988
c-------------------------------------------------------------------

      parameter (nmax=100,eps=5.e-6)
      parameter (zero=0.,one=1.,two=2.,three=3.,ten=10.)
      
      real sc,sb,ssi,sp,alk,t,s,ah
      real ph,co2,hco3,co3,pco2
      real g(4,3),sum(4)
      
      ta=t+273.15
      tac=ta/100.
      sal=s
      
      ak1=ten**(13.7201-0.031334*ta-3235.76/ta-1.3e-5*sal*ta
     .         +0.1032*sqrt(sal))
      a1=-5371.9645-1.671221*ta+128375.28/ta
      a2=2194.3055*log10(ta)-0.22913*sal-18.3802*log10(sal)
      a3=8.0944e-4*sal*ta+5617.11*log10(sal)/ta-2.136*sal/ta
      ak2=ten**(a1+a2+a3)

      akb=ten**(-9.26+0.00886*sal+0.01*t)

      aksi=4e-10

      akp1=2e-2
      akp2=exp(-9.039-1450/ta)
      akp3=exp(4.466-7276/ta)

      akw=exp(148.9802-13847.26/ta-23.6521*log(ta)
     .        -0.019813*sal+sqrt(sal)*(-79.2447+3298.72/ta
     .        +12.0408*log(ta)))
      fh=1.29-0.00204*ta+sal*sal*(4.61e-4-1.48e-6*ta)
      
      alfa=exp(-60.2409+9345.17/ta+23.3585*log(tac)+
     .         sal*(0.023517-0.023656*tac+0.0047036*tac*tac))
      
      
c     write(*,1000) ak1,ak2,akb,aksi
c1000  format(4(1pg20.12))
c      write(*,1000) akp1,akp2,akp3
c      write(*,1000) akw,fh
c      write(*,1000) alfa
      
      n=0
      a=one/ah
      sum(1)=sc
      sum(2)=sb
      sum(3)=ssi
      sum(4)=sp
      alka=alk
      
10    continue
      n=n+1
      al=a
      
      g(1,1)=one+ak1*a*(one+ak2*a)
      g(1,2)=ak1*(one+two*ak2*a)
      g(1,3)=two*ak1*ak2
      
      g(2,1)=one+akb*a
      g(2,2)=akb
      g(2,3)=zero
      
      g(3,1)=one+aksi*a
      g(3,2)=aksi
      g(3,3)=zero
      
      g(4,1)=one+akp1*a*(one+akp2*a*(one+akp3*a))
      g(4,2)=akp1*(one+akp2*a*(two+three*akp3*a))
      g(4,3)=akp1*akp2*(two+6*akp3*a)
      
      h=zero
      hs=zero
      do 20 l=1,4
c        write(*,1000) (g(l,j),j=1,3)
        h=h+sum(l)*g(l,2)/g(l,1)
        hs=hs+sum(l)*(g(l,3)*g(l,1)-g(l,2)*g(l,2))/(g(l,1)*g(l,1))
20    continue        

c      write(*,1000) h,hs
      hs=h+a*hs+1.e6*akw*fh+1.e6/(a*a*fh)
      h=a*(h+1.e6*akw*fh)-1.e6/(a*fh)-alka
      a=a-h/hs
c      write(9,'(i4,1pg20.12)') n,a
      if((abs((a-al)/al).gt.eps).and.(n.lt.nmax)) then
        goto 10
      endif
      
      if(n.ge.nmax) then
        write(*,*) 'No convergence in module CCHEMS, ah0=',ah
      endif
      
      ah=one/a
      hco3=sc/(1.+ak2/ah+ah/ak1)
      co3=hco3*ak2/ah
      co2=hco3*ah/ak1
      pco2=co2/alfa
      ph=-log10(ah/fh)

      return
      end
%}


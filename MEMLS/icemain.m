function [Tb,ev,eh]=icemain(ffile,tfile,par)
%   MEMLS
%   main program: brightness temperature of snowpack atop packice
%     ffile:  Inputfile with frequency array [GHz]
%     tfile:  Inputfile with incidence angle [deg]
%     ifile: Inputfile with layer- number, temp [K], vol.water content, density [kg/m3], thickness [cm], cor.length [mm],  salinity [per mil], ice/snow[1/0]
%     type:  ice first year = 3 or multiyear = 4
%   Version history:
%      1.0    wi 19.03.97
%      2.0    wi 25.09.97 adapted to the new subroutines
%      2.1    wi 15.10.97 adapted to the new subroutines
%      2.2    wi 04.11.97 adapted to bug in polmix
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland
%   Extension to sea ice by rtt dmi.dk
%   Revised by comparison to "Memls_mod.m" by tms space.dtu.dk

%   Constants:
%   sccho: type of scattering coefficient (11 recommended, 4=Born approximation spherical shape, 5=iborn for shells and spheres, now a function of type)
%   s0h/s0v : snow-ice reflectivity H/V-pol = 0
%   Tsky    : sky brightness temperature [K] = 0 (surface contribution to atmospheric total = radiometer sitting on ice)
%   Tgnd    : ground temperature [K] = -1.8 C (ocean beneath sea ice)
% pa col: [num,Ti,Wi,roi,di,pci,sal,si]
  Tsky=0;
  sccho=2;
  s0h=0;
  s0v=0;
  Tgnd=271.35;

  frq = [6.9,10.7,18.7,23.8,36.5]'; 
  tet = [55,55,55,55,55]';
  tet = (tet * pi) / 180;

  num = par(:,1);
  Ti = par(:,2);
  Wi = par(:,3);
  roi = par(:,4);
  di = par(:,5);
  pci = par(:,6);
  sal = par(:,7);
  si = par(:,8);
  type = par(:,9);
  format long

  N = length(num);
  
  
  if N == 0 
    return
  end
  roi = roi./1000;
  di = di./100;
  Tsky=0;
  %  loop over frequency
  
  for n=1:length(frq)
    freq=frq(n);
    teta=tet(n);

    % Absorption coefficient of all layers for calculating Teff
    [epsi,epsii] = ro2epsdd(roi,Ti,freq,2);

    %[epsi,epsii] = ro2epsd(roi,Ti,freq);
    [epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii);
    % Select dielectric scheme for FY or MY ice
    fy=(type==3);
    my=(type==4);
        
    [epsi,epsii] = sie(fy,sal,Ti,freq,epsi,epsii);
    
    [epsi,epsii] = mysie(my,roi,Ti,sal,freq,epsi,epsii);
   
    [gai] = abscoeff2(epsi,epsii,Ti,freq,Wi);

    ns = sqrt(epsi);
    tei = [asin(sin(teta)./ns);teta];
    
    dei = pfadi(tei,di);
    
    [sih,siv] = fresnelc(tei,[epsi;1]);
     
    [rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rtype,rsal] ...
        = slred2(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,type,sal);
    
    [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,sccho);
   
    %% All snow is new and old, recrystallized snow is ignored
    %% rtype = 1;
    rmeteo=(rtype==1);
    rrecry=(rtype==2);
    [gbih,gbiv,gs6,ga2i] = meteo_sc2(rroi,rTi,rpci,freq,rWi,rgai,rmeteo,gbih,gbiv,gs6,ga2i);
    [gbih,gbiv,gs6,ga2i] = recry_sc2(rroi,rTi,rpci,freq,rWi,rgai,rrecry,gbih,gbiv,gs6,ga2i);
    % Select scattering module for MY or FY ice
    %% Type of ice is a predefined scalar
    %% rtype = type;
    rfy=(rtype==3);
    rmy=(rtype==4);
     
    [gbih,gbiv,gs6,ga2i] = scice(rfy,gbih,gbiv,gs6,ga2i,rTi,rsal,freq, ...
                                 rpci);
   
    [gbih,gbiv,gs6,ga2i] = scice_my(rmy,gbih,gbiv,gs6,ga2i,rTi,rroi,freq,rpci,rsal);
    % Recompute 2-flux absorption coefficient (ga2i)
    
    [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,repsi,repsii,rroi,rTi,rpci,freq,rWi,rgai);

    [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6);
    
    rsih = [s0h;rsih];
    rsiv = [s0v;rsiv];
    [rsih,rsiv] = polmix(tscat,rsih,rsiv);
    rtei.*180./pi;
    
    [ri,ti] = rt(ga2i,gbih,rdei);

    Dh   = layer(ri,rsih,ti,rTi,Tgnd,Tsky);
    N = length(rroi);
    Tbh = (1-rsih(N+1))*Dh(N) + rsih(N+1)*Tsky
 
    [ri,ti]  = rt(ga2i,gbiv,rdei);
    Dv  = layer(ri,rsiv,ti,rTi,Tgnd,Tsky);
    Tbv = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*Tsky
    %Tbv
 
    Tb(n,:)=[Tbv,Tbh];
    
    eDh   = layer(ri,rsih,ti,rTi,Tgnd,0);
    eTbh = (1-rsih(N+1))*eDh(N);
    yh0(n) = eTbh;

    eDh = layer(ri,rsih,ti,rTi,Tgnd,100);
    eTbh = (1-rsih(N+1))*eDh(N) + rsih(N+1)*100;
    yh100(n) = eTbh;

    eDv   = layer(ri,rsiv,ti,rTi,Tgnd,0);
    eTbv = (1-rsiv(N+1))*eDv(N);
    yv0(n) = eTbv;

    eDv   = layer(ri,rsiv,ti,rTi,Tgnd,100);
    eTbv = (1-rsiv(N+1))*eDv(N) + rsiv(N+1)*100;
    yv100(n) = eTbv;

  end

  % calculate emissivities
  rv = (yv100 - yv0)./100;
  rh = (yh100 - yh0)./100;
  ev = 1 - rv;
  eh = 1 - rh;
end % end of function


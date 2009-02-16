%******************************
% Forward algorithm
% References: AMSR Ocean Algorithm; Frank J. Wentz, Thomas Meissner; Remote
%             Sensing Systems; Version 2; November 2; 2000.
% Dorthe Hofman-Bang
% 20th of September 2003
%*******************************
% Input: W,V,L,T_s,C_is,F_MY, 
%        Tb_MY,Tb_FY: Brightness temperature computed using  MEMLS  [6v,
%        6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h] [K]
%        T_memls: physical temperature of memls snow-ice model 
% Output: T_B[[6v, 6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h]]
%*******************************
function[T_B]=fw_fun3(p,ice_type)

  T_s       = p(1);
  V         = p(2);
  di_snow   = p(3);
  roi_snow  = p(4);
  pci_snow  = p(5);
  di_ice    = p(6);
  sal       = p(7);
 % FMY       = p8;

  snow_layernums=10;
  ice_type =1;

  if ice_type==1 %MY
    
    fid = fopen('input_MY', 'wt');
    MY_input = [100,T_s;0,roi_snow;di_ice,di_snow;0,pci_snow;sal,0];
    fprintf(fid,'%3.12f %3.12f %4.12f %4.15f %4.15f\n',MY_input);
    fclose(fid);

    discretize('input_MY',snow_layernums,'syntes.disc.MY',1);
    Tb_MY = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
                    'syntes.disc.MY',4);
   % emissivity_MY =Tb_MY/T_s;
    
   emissivity_MY = emi('syntes.disc.MY',4);
     emissivity_MY
  end

  if ice_type==2 %FY
    pci_FY = 0.8*3*Wi*0.333/4*4.14;
    Memls_FY = [T_sea, Wi_ice,roi_FY,di_ice,pci_FY,sal,1;T_s, Wi_snow,roi_snow,di_snow,pci_snow,0,0];
    fid = fopen('input_FY', 'wt');
    fprintf(fid,'%3.12f %3.12f %4.12f %4.15f %4.15f %3.15f %2d\n', Memls_FY');
    fclose(fid);

    discretize('input_FY',snow_layernums,'syntes.disc.FY');
    Tb_FY = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
                    'syntes.disc.FY',3);
    emissivity_FY = emi('syntes.disc.FY',3);
  end

  theta_d=55;
  theta_r=55*3.1415926535897/180;      % Incidence angle [rad] NB!!  

  %************************
  % Constants
  % Columns:[[6.93],  [10.65],    [18.70],   [23.80],    [36.50],   [50.30],    [52.80],    [89.00]] (GHz)
  %************************
  Freq=[ 6.93;
         10.65;
         18.70;
         23.80;
         36.50]; %AMSR frequencies

  %Model Coefficients for the atmosphere
              %b0,        b1,        b2,        b3,        b4,       b5,       b6,       b7,      a01,     a02,      aV1,      aV2;
  Mc_atm=[  239.50E+0, 213.92E-2, -460.60E-4, 457.11E-6, -16.84E-7, 0.50E+0, -0.11E+0, -0.21E-2, 8.34E-3, -0.48E-4, 0.07E-3, 0.00E-5;
            239.51E+0, 225.19E-2, -446.86E-4, 391.82E-6, -12.20E-7, 0.54E+0, -0.12E+0, -0.34E-2, 9.08E-3, -0.47E-4, 0.18E-3, 0.00E-5;
            240.24E+0, 298.88E-2, -725.93E-6, 814.50E-6, -36.07E-7, 0.61E+0, -0.16E+0, -1.69E-2, 12.15E-3, -0.61E-4, 1.73E-3, -0.05E-5;
            241.69E+0, 310.32E-2, -814.29E-4, 998.93E-6, -48.37E-7, 0.20E+0, -0.20E+0, -5.21E-2, 15.75E-3, -0.87E-4, 5.14E-3, 0.19E-5;
            239.45E+0, 254.41E-2, -512.84E-4, 452.02E-6, -14.36E-7, 0.58E-0, -0.57E+0, -2.38E-2, 40.06E-3, -2.00E-4, 1.88E-3, 0.09E-5]; 

  % Coefficients for Rayleigh Absorption and Mie Scattering
            %aL1,       aL2,         aL3,      aL4,       aL5;
  Mc_abs=[0.0078,    0.0303,    0.0007,         0,    1.2216;
          0.0183,    0.0298,    0.0027,    0.0060,    1.1795;
          0.0556,    0.0288,    0.0113,    0.0040,    1.0636;
          0.0891,    0.0281,    0.0188,    0.0020,    1.0220;
          0.2027,    0.0261,    0.0425,   -0.0020,    0.9546];

  % Model coefficients for Geometric optics
           %v-pol r0,  h-pol r0, v-pol r1,    h-pol r1,  v-pol r2,   h-pol r2,      v-pol r3,    h-pol r3;
  Mc_geo =[-0.00027,   0.00054,  -0.000021,   0.000032,  -0.000021,  -0.00002526,         0,               0;
           -0.00032,   0.00072,  -0.000029,   0.000044,  -0.000021,  -0.00002894,   0.00000008,  -0.00000002;
           -0.00049,   0.00113,  -0.000053,   0.000070,  -0.000021,  -0.00003690,   0.00000031,  -0.00000012;
           -0.00063,   0.00139,  -0.000070,   0.000085,  -0.000021,  -0.00004195,   0.00000041,  -0.00000020;
           -0.00101,   0.00191,  -0.000105,   0.000112,  -0.000021,  -0.00005451,   0.00000045,  -0.00000036];

  % Model coefficients m
       %v-pol m1,  h-pol m1,  v-pol m2, h-pol m2;
  Mc_m =[0.00020,   0.00200,   0.00690,   0.00600;
         0.00020,   0.00200,   0.00690,   0.00600;
         0.00140,   0.00293,   0.00736,   0.00656;
         0.00178,   0.00308,   0.00730,   0.00660;
         0.00257,   0.00329,   0.00701,   0.00660];

  T_C=2.7;     %Background radiation [K]

%   T_BV_FY= Tb_FY(:,1);
%   T_BH_FY= Tb_FY(:,2);
  T_BV_MY = Tb_MY(:,1);
  T_BH_MY = Tb_MY(:,2);

  T_S_mix=T_s;    %Temperature of a mixed surface (ice and open water)

%  T_L=(T_S_mix+273)/2;%Temperature of the water droplets [Kelvin]; approximation: mean temp of surface and freezing level 

  %************************
  % Model
  %************************
  % Model for the Atmosphere

  if V<=48
    T_V=273.16+0.8337*V-3.029E-5*V^3.33; %equation (27a)
  else
    T_V=301.16;                          %equation (27b)
  end

  temp_test=abs(T_S_mix-T_V);
  if temp_test<=20
    sig_TS_TV=1.05*(T_S_mix-T_V)*(1-((T_S_mix-T_V)^2)/1200); %equation (27c)
  else
    sig_TS_TV=sign(T_S_mix-T_V)*14;                      %equation (27d)
  end

  T_D=Mc_atm(:,1)+Mc_atm(:,2)*V+Mc_atm(:,3)*V^2+Mc_atm(:,4)*V^3+Mc_atm(:,5)*V^4+Mc_atm(:,6)*sig_TS_TV;    %equation (26a)

  T_U=T_D+Mc_atm(:,7)+Mc_atm(:,8)*V;     %equation (26b)

  A_0=Mc_atm(:,9)+Mc_atm(:,10).*(T_D-270);        %equation (28)
  A_V=Mc_atm(:,11)*V+Mc_atm(:,12)*V^2;            %equation (29)

  % A_L=Mc_abs(:,1).*(1-Mc_abs(:,2)*(T_L-283))*L;   %equation (33) (no rain)

  tau=exp((-1/cos(theta_r))*(A_0+A_V));    %Total transmittance; equation (22)

  T_BU=T_U.*(1-tau);   %The upwelling effective air temperature; equation (24a)
  T_BD=T_D.*(1-tau);   %The downwelling effective air temperature; equation (24b)


  %Emissivity for mixed surface
%  
%   E_eff_H=C_FY*E_FY(1:5,2)+C_MY*E_MY(1:5,2);
%    E_eff_V=C_FY*E_FY(1:5,1)+C_MY*E_MY(1:5,1);
    
      E_eff_H=emissivity_MY(1:5,2);
      E_eff_V=emissivity_MY(1:5,1);

  R_eff_H=1-E_eff_H;   %Reflection coefficient for mixed surface
  R_eff_V=1-E_eff_V;


  %************************
  % Output
  %************************
%   T_BS_H= C_MY*Tb_MY(1:5,2) + C_FY*Tb_FY(1:5,2);
%   T_BS_V= C_MY*Tb_MY(1:5,1) + C_FY*Tb_FY(1:5,1);
 T_BS_H= Tb_MY(1:5,2) ;
  T_BS_V= Tb_MY(1:5,1);

  T_BOmegaH = (T_BD+T_C*tau).*R_eff_H;
  T_BOmegaV = (T_BD+T_C*tau).*R_eff_V;

  T_BH=(T_BU+(tau.*(T_BS_H+T_BOmegaH)))';     %Equation (10); Horizontal pol.
  T_BV=(T_BU+(tau.*(T_BS_V+T_BOmegaV)))';     %Equation (10); Vertical pol.


  T_B=[T_BV(1); T_BH(1); T_BV(2); T_BH(2); T_BV(3); T_BH(3); T_BV(4); T_BH(4); ...
       T_BV(5); T_BH(5)]';


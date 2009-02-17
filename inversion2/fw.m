%******************************
% Forward algorithm
%%*******************************
% Input: 
% Output: T_B[[6v, 6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h]]
%*******************************
function[T_B]=fw(p)
V=2;
  theta_d=55;
theta_r=55*3.1415926535897/180;      % Incidence angle [rad] NB!!  
%************************
% Constants
% Columns:[[6.93],  [10.65],    [18.70],   [23.80],    [36.50],   [50.30],    [52.80],    [89.00]] (GHz)
%************************
T_snow=p(1,1);
di_snow=p(2,1);
pci_snow=p(3,1);
roi_snow=p(4,1);
C_FY=p(5,1);


%[T_snow_FY,T_snow_MY,T_ice_FY,T_ice_MY]  = phytemp(p,1);
T_ice_FY = 263;
T_ice_MY = 255.5;
FY= [1,  T_ice_FY,  0.1,  920     ,      100,     0.17,  7,  1, 3;
     2,  T_snow,  0,  roi_snow,  di_snow, pci_snow,  0 , 0, 1];
MY= [1,  T_ice_MY,  0.15,  900     ,      300,     0.25, 2  1, 4;
     2,  T_snow,  0,  roi_snow,  di_snow, pci_snow,  0 , 0, 2];


Freq=[ 6.93;
      10.65;
      18.70;
      23.80;
      36.50];%, 50.30, 52.80, 89.00]'; %AMSR frequencies
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

T_C=2.7;     %Background radiation [K
%%%%%%%%%%%%%%%%%%% Bias 

 Tb_FY_bias= [0.48,-5.86;3.96,-5.61;-1.37,-8.77;-2.58, -10.69;-3.99, -8.35];
 Tb_MY_bias=[2.49,0.85;6.49,2.87;5.93,2.74;5.22,2.49;2.48,3.79];
  
  
if C_FY<0
     C_MY=0;
end
if C_FY>1
     C_FY=1;
end
C_MY = 1-C_FY;
cd ../MEMLS
[Tb_FY,ev_FY,eh_FY]=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);
Tb_FY = Tb_FY-Tb_FY_bias;
[Tb_MY,ev_MY,eh_MY]=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY);
cd ../inversion2

Tb_V_FY=Tb_FY(:,1);
Tb_H_FY=Tb_FY(:,2);
Tb_V_MY=Tb_MY(:,1);
Tb_H_MY=Tb_MY(:,2);

T_sea=271.35;      %[K] 
T_S_FY = 2*T_snow+T_sea-2*T_ice_FY;
T_S_MY = 2*T_snow+T_sea-2*T_ice_MY;
T_S_mix = C_FY*T_S_FY + C_MY*T_S_MY; 

%T_S_mix = T_snow;
T_L=(T_S_mix+273)/2;%Temperature of the water droplets [Kelvin];
                        %approximation: mean temp of surface and freezing
                        %level 

%************************
% Model
%************************
% Model for the Atmosphere

if V<=48
    T_V=273.16+0.8337*V-3.029E-5*V^3.33; %equation (27a)
else
    T_V=301.16;                          %equation (27b)
end
%%%%%%%%%%%%%%%%%%%%% V
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
%%%%%%%%%%%%%%%%%%%% H
%Emissivity for mixed surface
E_eff_H=C_FY*eh_FY+C_MY*eh_MY;
E_eff_V=C_FY*ev_FY+C_MY*ev_MY;

R_eff_H=1-E_eff_H;   %Reflection coefficient for mixed surface
R_eff_V=1-E_eff_V;

%************************
% Output
%************************
T_BS_H= C_MY*Tb_H_MY+C_FY*Tb_H_FY;
T_BS_V= C_MY*Tb_V_MY+C_FY*Tb_V_FY;

T_BOmegaH = (T_BD+T_C*tau).*R_eff_H';
T_BOmegaV = (T_BD+T_C*tau).*R_eff_V';
T_BV=(T_BU+(tau.*(T_BS_V+T_BOmegaV)))';     %Equation (11); Vertical pol.
T_BH=(T_BU+(tau.*(T_BS_H+T_BOmegaH)))';     %Equation (10); Horizontal pol.

%************************

T_B=[T_BV(1); T_BH(1); T_BV(2); T_BH(2); T_BV(3); T_BH(3); T_BV(4); T_BH(4); ...
     T_BV(5); T_BH(5)]';


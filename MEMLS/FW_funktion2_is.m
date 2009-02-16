%******************************
% Forward algorithm
% References: AMSR Ocean Algorithm; Frank J. Wentz, Thomas Meissner; Remote
%             Sensing Systems; Version 2; November 2; 2000.
% Dorthe Hofman-Bang
% 20th of September 2003
%*******************************
% Input: W,V,L,T_s,C_is,F_MY 
% Output: T_B[[6v, 6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h]]
%*******************************
function[T_B]=FW_funktion2_is(W, V, L, T_ow, T_is, C_is,F_MY)

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
      36.50];%, 50.30, 52.80, 89.00]'; %AMSR frequencies

%Model Coefficients for the atmosphere
%Mc_atm=[  239.50E+0,  239.51E+0,  240.24E+0,  241.69E+0,  239.45E+0; %,  242.10E+0,  245.87E+0,  242.58E+0; %b0 (K)
%          213.92E-2,  225.19E-2,  298.88E-2,  310.32E-2,  254.41E-2; %,  229.17E-2,  250.61E-2,  302.33E-2; %b1 (K mm?1)
%         -460.60E-4, -446.86E-4, -725.93E-4, -814.29E-4, -512.84E-4; %, -508.05E-4, -627.89E-4, -749.76E-4; %b2 (K mm-2)
%          457.11E-6,  391.82E-6,  814.50E-6,  998.93E-6,  452.02E-6; %,  536.90E-6,  759.62E-6,  880.66E-6; %b3 (K mm-3)
%          -16.84E-7,  -12.20E-7,  -36.07E-7,  -48.37E-7,  -14.36E-7; %,  -22.07E-7,  -36.06E-7,  -40.88E-7; %b4 (K mm-4)
%            0.50E+0,    0.54E+0,    0.61E+0,    0.20E+0,    0.58E+0; %,    0.52E+0,    0.53E+0,    0.62E+0; %b5
%           -0.11E+0,   -0.12E+0,   -0.16E+0,   -0.20E+0,   -0.57E+0; %,   -4.59E+0,  -12.52E+0,   -0.57E+0; %b6 (K)
%           -0.21E-2,   -0.34E-2,   -1.69E-2,   -5.21E-2,   -2.38E-2; %,   -8.78E-2,  -23.26E-2,   -8.07E-2; %b7 (K mm-1)
%            8.34E-3,    9.08E-3,   12.15E-3,   15.75E-3,   40.06E-3; %,  353.72E-3, 1131.76E-3,   53.35E-3; %aO1
%           -0.48E-4,   -0.47E-4,   -0.61E-4,   -0.87E-4,   -2.00E-4; %,  -13.79E-4,   -2.26E-4,   -1.18E-4; %aO2 (K-1)
%            0.07E-3,    0.18E-3,    1.73E-3,    5.14E-3,    1.88E-3; %,    2.91E-3,    3.17E-3,    8.78E-3; %aV1 (mm-1)
%            0.00E-5,    0.00E-5,   -0.05E-5,    0.19E-5,    0.09E-5]'; %,    0.24E-5,    0.27E-5,    0.80E-5]'; %aV2 (mm-2)

            %b0,        b1,        b2,        b3,        b4,       b5,       b6,       b7,      a01,     a02,      aV1,      aV2;
Mc_atm=[  239.50E+0, 213.92E-2, -460.60E-4, 457.11E-6, -16.84E-7, 0.50E+0, -0.11E+0, -0.21E-2, 8.34E-3, -0.48E-4, 0.07E-3, 0.00E-5;
          239.51E+0, 225.19E-2, -446.86E-4, 391.82E-6, -12.20E-7, 0.54E+0, -0.12E+0, -0.34E-2, 9.08E-3, -0.47E-4, 0.18E-3, 0.00E-5;
          240.24E+0, 298.88E-2, -725.93E-6, 814.50E-6, -36.07E-7, 0.61E+0, -0.16E+0, -1.69E-2, 12.15E-3, -0.61E-4, 1.73E-3, -0.05E-5;
          241.69E+0, 310.32E-2, -814.29E-4, 998.93E-6, -48.37E-7, 0.20E+0, -0.20E+0, -5.21E-2, 15.75E-3, -0.87E-4, 5.14E-3, 0.19E-5;
          239.45E+0, 254.41E-2, -512.84E-4, 452.02E-6, -14.36E-7, 0.58E-0, -0.57E+0, -2.38E-2, 40.06E-3, -2.00E-4, 1.88E-3, 0.09E-5]; 

% Coefficients for Rayleigh Absorption and Mie Scattering
%Mc_abs=[0.0078, 0.0183, 0.0556, 0.0891,  0.2027; %,  0.3682,  0.4021,  0.9693; %aL1
%        0.0303, 0.0298, 0.0288, 0.0281,  0.0261; %,  0.0236,  0.0231,  0.0146; %aL2
%        0.0007, 0.0027, 0.0113, 0.0188,  0.0425; %,  0.0731,  0.0786,  0.1506; %aL3
%        0.0000, 0.0060, 0.0040, 0.0020, -0.0020; %, -0.0020, -0.0020, -0.0020; %aL4
%        1.2216, 1.1795, 1.0636, 1.0220,  0.9546]'; %,  0.8983,  0.8943,  0.7961]'; %aL5

        %aL1,       aL2,         aL3,      aL4,       aL5;
Mc_abs=[0.0078,    0.0303,    0.0007,         0,    1.2216;
        0.0183,    0.0298,    0.0027,    0.0060,    1.1795;
        0.0556,    0.0288,    0.0113,    0.0040,    1.0636;
        0.0891,    0.0281,    0.0188,    0.0020,    1.0220;
        0.2027,    0.0261,    0.0425,   -0.0020,    0.9546];

% Model coefficients for Geometric optics
%Mc_geo=[ -0.27E-3,  -0.32E-3,  -0.49E-3,  -0.63E-3,  -1.01E-3; %, -1.20E-3, -1.23E-03, -1.53E-3; %v-pol r0
%          0.54E-3,   0.72E-3,   1.13E-3,   1.39E-3,   1.91E-3; %,  1.97E-3,  1.97E-03,  2.02E-3; %h-pol r0
%        -0.21E-4,  -0.29E-4,  -0.53E-4,  -0.70E-4,  -1.05E-4; %, -1.12E-4, -1.13E-04, -1.16E-4; %v-pol r1
%          0.32E-4,   0.44E-4,   0.70E-4,   0.85E-4,   1.12E-4; %,  1.18E-4,  1.19E-04,  1.30E-4; %h-pol r1
%         -2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5; %, -2.10E-5, -2.10E-05, -2.10E-5; %v-pol r2 empirically
%        -25.26E-6, -28.94E-6, -36.90E-6, -41.95E-6, -54.51E-6; %, -5.50E-5, -5.50E-5,  -5.50E-5; %h-pol r2 empirically
%          0.00E-6,   0.08E-6,   0.31E-6,   0.41E-6,   0.45E-6; %,  0.35E-6,  0.32E-06, -0.09E-6; %v-pol r3
%          0.00E-6,  -0.02E-6,  -0.12E-6,  -0.20E-6,  -0.36E-6]'; %, -0.43E-6, -0.44E-06, -0.46E-6]'; %h-pol r3

         %v-pol r0,  h-pol r0, v-pol r1,    h-pol r1,  v-pol r2,   h-pol r2,      v-pol r3,    h-pol r3;
Mc_geo =[-0.00027,   0.00054,  -0.000021,   0.000032,  -0.000021,  -0.00002526,         0,               0;
         -0.00032,   0.00072,  -0.000029,   0.000044,  -0.000021,  -0.00002894,   0.00000008,  -0.00000002;
         -0.00049,   0.00113,  -0.000053,   0.000070,  -0.000021,  -0.00003690,   0.00000031,  -0.00000012;
         -0.00063,   0.00139,  -0.000070,   0.000085,  -0.000021,  -0.00004195,   0.00000041,  -0.00000020;
         -0.00101,   0.00191,  -0.000105,   0.000112,  -0.000021,  -0.00005451,   0.00000045,  -0.00000036];

% Model coefficients m
%Mc_m=[  0.00020, 0.00020, 0.00140, 0.00178, 0.00257; %, 0.00260, 0.00260, 0.00260; %v-pol m1
%        0.00200, 0.00200, 0.00293, 0.00308, 0.00329; %, 0.00330, 0.00330, 0.00330; %h-pol m1
%        0.00690, 0.00690, 0.00736, 0.00730, 0.00701; %, 0.00700, 0.00700, 0.00700; %v-pol m2
%        0.00600, 0.00600, 0.00656, 0.00660, 0.00660]'; %, 0.00660, 0.00660, 0.00660]'; %h-pol m2

      %v-pol m1,  h-pol m1,  v-pol m2, h-pol m2;
Mc_m =[0.00020,   0.00200,   0.00690,   0.00600;
       0.00020,   0.00200,   0.00690,   0.00600;
       0.00140,   0.00293,   0.00736,   0.00656;
       0.00178,   0.00308,   0.00730,   0.00660;
       0.00257,   0.00329,   0.00701,   0.00660];

T_C=2.7;     %Background radiation [K]

Emissivitet_FY_V=[0.9204; 0.9127; 0.9373; 0.9409; 0.9347]; % 0; 0; 0;];
Emissivitet_FY_H=[0.7502; 0.7738; 0.8314; 0.8490; 0.8600]; % 0; 0; 0;];
Emissivitet_MY_V=[0.9692; 0.9284; 0.8843; 0.8554; 0.7813]; % 0; 0; 0;];
Emissivitet_MY_H=[0.8651; 0.8356; 0.7917; 0.7792; 0.7248]; % 0; 0; 0;];

%****************************************
%Information om is
%****************************************
% Columns:[[6.93],  [10.65],    [18.70],   [23.80],    [36.50],   [50.30],    [52.80],    [89.00]] (GHz)
               
%Emissivitet_FY_V=[0.9204; 0.9127; 0.9373; 0.9409; 0.9347; 0; 0; 0;];
%Emissivitet_FY_H=[0.7502; 0.7738; 0.8314; 0.8490; 0.8600; 0; 0; 0;];
%Emissivitet_MY_V=[0.9692; 0.9284; 0.8843; 0.8554; 0.7813; 0; 0; 0;];
%Emissivitet_MY_H=[0.8651; 0.8356; 0.7917; 0.7792; 0.7248; 0; 0; 0;];
         
% Emissiviteter fra Bog
% Emissivitet_FY_H=[0.84;	0.876;	0.888;	0.91;	0.913; 0; 0; 0];
% Emissivitet_FY_V=[0.9;	0.924;	0.941;	0.96;	0.955; 0; 0; 0];
% Emissivitet_MY_H=[0.82;	0.817;	0.78;	0.635;	0.706; 0; 0; 0];
% Emissivitet_MY_V=[0.925;	0.89;	0.85;	0.787;	0.764; 0; 0; 0];


%C_is=C_MY+C_FY;  %Total iskoncentration
if F_MY<0
     F_MY=0;
end
if F_MY>1
     F_MY=1;
end
% if C_is>1
%     C_is=1;
% end
% if C_is<0
%     C_is=0;
% end


C_MY=C_is*F_MY; %Multi Year ice concentration
C_FY=C_is-C_MY; %First Year ice concentration
C_ow=1-C_is;    %Open water concentration

%s\ufffdt overflade temperatur
%if T_ow>273.16
    %T_is=273.16;
    %else
   % T_is=0.4*T_ow+0.6*272;
    %end

%if C_is>0.05
% T_ow=273.16;
%else
    T_ow=T_ow;
    %end

T_BV_FY=T_is*Emissivitet_FY_V;  %Brigthness temperature from FY is (vertical pol)
T_BH_FY=T_is*Emissivitet_FY_H;  %Brigthness temperature from FY is (horizontal pol)

T_BV_MY=T_is*Emissivitet_MY_V;  %Brigthness temperature from MY is (vertical pol)
T_BH_MY=T_is*Emissivitet_MY_H;  %Brigthness temperature from MY is (horizontal pol)

T_S_mix=C_is*T_is+C_ow*T_ow;    %Temperature of a mixed surface (ice and open water)

T_L=(T_S_mix+273)/2;%Temperature of the water droplets [Kelvin]; approximation: mean temp of surface and freezing level 

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

A_L=Mc_abs(:,1).*(1-Mc_abs(:,2)*(T_L-283))*L;   %equation (33) (no rain)

tau=exp((-1/cos(theta_r))*(A_0+A_V+A_L));    %Total transmittance; equation (22)

T_BU=T_U.*(1-tau);   %The upwelling effective air temperature; equation (24a)
T_BD=T_D.*(1-tau);   %The downwelling effective air temperature; equation (24b)

%********************************
%Dielectric Constant of Sea-water 
t_ow=T_ow-273.16;      % Surface temperature T_ow [deg Celcius]
epsilon_R=4.44;      % Dielectric constant at inf. freq.
s=35;                % Salinity in parts per thousand; Typical values 30-33ppt
ny=0.012;            % Spread factor
light_speed=3.00E10; % Speed of light, [cm/s]
lambda=(light_speed./(Freq*1E9));    % Wavelength, [cm]
%free_space_permittivity=8.854E-12

epsilon_S=(87.90*exp(-0.004585*t_ow))*(exp(-3.45E-3*s+4.69E-6*s^2+1.36E-5*s*t_ow));       %equation (36) og (43)
lambda_R=(3.30*exp(-0.0346*t_ow+0.00017*t_ow^2))-(6.54E-3*(1-3.06E-2*t_ow+2.0E-4*t_ow^2)*s);   %equation (38) og (44)

C=0.5536*s;      %equation (41)
delta_t=25-t_ow;  %equation (42)
qsi=2.03E-2+1.27E-4*delta_t+2.46E-6*delta_t^2-C*(3.34E-5-4.60E-7*delta_t+4.60E-8*delta_t^2); %equation (40)
sigma=3.39E9*C^0.892*exp(-delta_t*qsi);      %equation (39)

epsilon=epsilon_R+((epsilon_S-epsilon_R)./(1+((j*lambda_R)./lambda).^(1-ny)))-((2*j*sigma.*lambda)./light_speed);   %equation (35)

rho_H=(cos(theta_r)-sqrt(epsilon-sin(theta_r)^2))./(cos(theta_r)+sqrt(epsilon-sin(theta_r)^2));   %equation (45b); horizontal pol.
rho_V=(epsilon*cos(theta_r)-sqrt(epsilon-sin(theta_r)^2))./(epsilon*cos(theta_r)+sqrt(epsilon-sin(theta_r)^2));   %equation (45a); vertical pol.

R_0H=abs(rho_H).^2;   %equation (46); horizontal pol.
R_0V=abs(rho_V).^2+(4.887E-8-6.108E-8*(T_ow-273)^3);   %equation (46)+correction vertical pol.

%********************************
% The wind-roughened Sea Surface

R_geoH=R_0H-(Mc_geo(:,2)+Mc_geo(:,4)*(theta_d-53)+Mc_geo(:,6)*(T_ow-288)+Mc_geo(:,8)*(theta_d-53)*(T_ow-288))*W; %equation(57); horizontal pol.
R_geoV=R_0V-(Mc_geo(:,1)+Mc_geo(:,3)*(theta_d-53)+Mc_geo(:,5)*(T_ow-288)+Mc_geo(:,7)*(theta_d-53)*(T_ow-288))*W; %equation(57); vertical pol.

%horizontal polarisation
W_1=7;
W_2=12;
if W<W_1
    F_H=Mc_m(:,2)*W;       %equation (60a)
elseif (W>=W_1 & W<=W_2)
    F_H=Mc_m(:,2)*W+0.5*(Mc_m(:,4)-Mc_m(:,2))*(W-W_1)^2/(W_2-W_1);   %equation (60b)
else
    F_H=Mc_m(:,4)*W-0.5*(Mc_m(:,4)-Mc_m(:,2))*(W_2+W_1); %equation (60c)
end

%vertikal polarisation
W_1=3;
W_2=12;
if W<W_1
    F_V=Mc_m(:,1)*W;   %equation (60a)
elseif (W>=W_1 & W<=W_2)
    F_V=Mc_m(:,1)*W+0.5*(Mc_m(:,3)-Mc_m(:,1))*(W-W_1)^2/(W_2-W_1);   %equation (60b)
else
    F_V=Mc_m(:,3)*W-0.5*(Mc_m(:,3)-Mc_m(:,1))*(W_2+W_1); %equation (60c)
end

R_H=(1-F_H).*R_geoH;   %equation (49) Composite reflektivity for open water; horisontal pol
R_V=(1-F_V).*R_geoV;   %equation (49) Composite reflektivity for open water; vertical pol

E_H=1-R_H;       % surface emissivity for open water; Horizontal pol.; Equation (8); 
E_V=1-R_V;       % Surface emissivity for open water; Vertical pol; Equation (8); 

%Emissivity for mixed surface
E_eff_H=C_ow*E_H+C_FY*Emissivitet_FY_H+C_MY*Emissivitet_MY_H;
E_eff_V=C_ow*E_V+C_FY*Emissivitet_FY_V+C_MY*Emissivitet_MY_V;

R_eff_H=1-E_eff_H;   %Reflection coefficient for mixed surface
R_eff_V=1-E_eff_V;

%********************************
% Atmospheric Radiation Scattered by the Sea Surface
% Calculatet Delta_S2
Delta_S2(5,1)=5.22E-3*W;
Delta_S2(1:4,1)=5.22E-3*(1-0.00748*(37-Freq(1:4)).^1.3)*W;

for n=1:5
    if Delta_S2(n)>0.069        %equation (56a)
        Delta_S2(n)=0.069;      %equation (56b)
    end
end

term=Delta_S2-70*Delta_S2.^3; %Term for equation (62a+b)

OmegaH=(6.2-0.001*(37-Freq).^2).*term.*tau.^2.0;  %equation (62a); horizontal pol.
OmegaV=(2.5+0.018*(37-Freq)).*term.*tau.^3.4;     %equation (62b); vertical pol.

T_BOmegaH=((1+OmegaH).*(1-tau).*(T_D-T_C)+T_C).*R_eff_H;    %equation (61) horizontal pol.
T_BOmegaV=((1+OmegaV).*(1-tau).*(T_D-T_C)+T_C).*R_eff_V;    %equation (61) vertical pol.

%************************
% Output
%************************
T_BH_ow=E_H*T_ow;    %Horizontal brightnesstemperature from open water
T_BV_ow=E_V*T_ow;    %Vertikal brightnesstemperature from opent water

T_BH_overflade=C_ow*T_BH_ow+C_FY*T_BH_FY+C_MY*T_BH_MY;  %Horizontal brigthnesstemperature from mixed surface
T_BV_overflade=C_ow*T_BV_ow+C_FY*T_BV_FY+C_MY*T_BV_MY;  %Vertikal brigthnesstemperature from mixed surface

T_BH=(T_BU+(tau.*(T_BH_overflade+T_BOmegaH)))';     %Equation (10); Horizontal pol.
T_BV=(T_BU+(tau.*(T_BV_overflade+T_BOmegaV)))';     %Equation (10); Vertical pol.

%************************
% sort output (exklude 53GHz)

% Alle frekvenser [6v, 6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h, 50v, 50h ,53v, 53h, 89v, 89h]
%T_B=[T_BV(1), T_BH(1), T_BV(2), T_BH(2), T_BV(3), T_BH(3), T_BV(4), T_BH(4), T_BV(5), T_BH(5), T_BV(6), T_BH(6), T_BV(7), T_BH(7), T_BV(8), T_BH(8)]

% frekvenser i Amerikansk AMSR [6v, 6h, 10v, 10h, 18v, 18h, 23v, 23h, 37v, 37h]
T_B=[T_BV(1); T_BH(1); T_BV(2); T_BH(2); T_BV(3); T_BH(3); T_BV(4); T_BH(4); T_BV(5); T_BH(5)];

function test_ice_orbit1_funktion
%global Freq;
%global Mc_atm;
%global Mc_abs;
%global Mc_geo;
%global Mc_m;
%global T_C;
%global Emissivitet_FY_V;
%global Emissivitet_FY_H;
%global Emissivitet_MY_V;
%global Emissivitet_MY_H;

%Freq=[6.93, 10.65, 18.70, 23.80, 36.50]';%, 50.30, 52.80, 89.00]'; %AMSR frequencies

%Model Coefficients for the atmosphere
%Mc_atm=[  239.50E+0,  239.51E+0,  240.24E+0,  241.69E+0,  239.45E+0; %,  242.10E+0,  245.87E+0,  242.58E+0; %b0 (K)
%	  213.92E-2,  225.19E-2,  298.88E-2,  310.32E-2,  254.41E-2; %,  229.17E-2,  250.61E-2,  302.33E-2; %b1 (K mm?1)
%	 -460.60E-4, -446.86E-4, -725.93E-4, -814.29E-4, -512.84E-4; %, -508.05E-4, -627.89E-4, -749.76E-4; %b2 (K mm-2)
%	  457.11E-6,  391.82E-6,  814.50E-6,  998.93E-6,  452.02E-6; %,  536.90E-6,  759.62E-6,  880.66E-6; %b3 (K mm-3)
%	  -16.84E-7,  -12.20E-7,  -36.07E-7,  -48.37E-7,  -14.36E-7; %,  -22.07E-7,  -36.06E-7,  -40.88E-7; %b4 (K mm-4)
%	    0.50E+0,    0.54E+0,    0.61E+0,    0.20E+0,    0.58E+0; %,    0.52E+0,    0.53E+0,    0.62E+0; %b5
%	   -0.11E+0,   -0.12E+0,   -0.16E+0,   -0.20E+0,   -0.57E+0; %,   -4.59E+0,  -12.52E+0,   -0.57E+0; %b6 (K)
%	   -0.21E-2,   -0.34E-2,   -1.69E-2,   -5.21E-2,   -2.38E-2; %,   -8.78E-2,  -23.26E-2,   -8.07E-2; %b7 (K mm-1)
%	    8.34E-3,    9.08E-3,   12.15E-3,   15.75E-3,   40.06E-3; %,  353.72E-3, 1131.76E-3,   53.35E-3; %aO1
%	   -0.48E-4,   -0.47E-4,   -0.61E-4,   -0.87E-4,   -2.00E-4; %,  -13.79E-4,   -2.26E-4,   -1.18E-4; %aO2 (K-1)
%	    0.07E-3,    0.18E-3,    1.73E-3,    5.14E-3,    1.88E-3; %,    2.91E-3,    3.17E-3,    8.78E-3; %aV1 (mm-1)
%	    0.00E-5,    0.00E-5,   -0.05E-5,    0.19E-5,    0.09E-5]'; %,    0.24E-5,    0.27E-5,    0.80E-5]'; %aV2 (mm-2)

% Coefficients for Rayleigh Absorption and Mie Scattering
%Mc_abs=[0.0078, 0.0183, 0.0556, 0.0891,  0.2027; %,  0.3682,  0.4021,  0.9693; %aL1
%        0.0303, 0.0298, 0.0288, 0.0281,  0.0261; %,  0.0236,  0.0231,  0.0146; %aL2
%	0.0007, 0.0027, 0.0113, 0.0188,  0.0425; %,  0.0731,  0.0786,  0.1506; %aL3
%	0.0000, 0.0060, 0.0040, 0.0020, -0.0020; %, -0.0020, -0.0020, -0.0020; %aL4
%	1.2216, 1.1795, 1.0636, 1.0220,  0.9546]'; %,  0.8983,  0.8943,  0.7961]'; %aL5

% Model coefficients for Geometric optics
%Mc_geo=[ -0.27E-3,  -0.32E-3,  -0.49E-3,  -0.63E-3,  -1.01E-3; %, -1.20E-3, -1.23E-03, -1.53E-3; %v-pol r0
%	  0.54E-3,   0.72E-3,   1.13E-3,   1.39E-3,   1.91E-3; %,  1.97E-3,  1.97E-03,  2.02E-3; %h-pol r0
%	 -0.21E-4,  -0.29E-4,  -0.53E-4,  -0.70E-4,  -1.05E-4; %, -1.12E-4, -1.13E-04, -1.16E-4; %v-pol r1
%	  0.32E-4,   0.44E-4,   0.70E-4,   0.85E-4,   1.12E-4; %,  1.18E-4,  1.19E-04,  1.30E-4; %h-pol r1
%	 -2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5; %, -2.10E-5, -2.10E-05, -2.10E-5; %v-pol r2 empirically
%        -25.26E-6, -28.94E-6, -36.90E-6, -41.95E-6, -54.51E-6; %, -5.50E-5, -5.50E-5,  -5.50E-5; %h-pol r2 empirically
%	  0.00E-6,   0.08E-6,   0.31E-6,   0.41E-6,   0.45E-6; %,  0.35E-6,  0.32E-06, -0.09E-6; %v-pol r3
%	  0.00E-6,  -0.02E-6,  -0.12E-6,  -0.20E-6,  -0.36E-6]'; %, -0.43E-6, -0.44E-06, -0.46E-6]'; %h-pol r3

% Model coefficients m
%Mc_m=[  0.00020, 0.00020, 0.00140, 0.00178, 0.00257; %, 0.00260, 0.00260, 0.00260; %v-pol m1
%        0.00200, 0.00200, 0.00293, 0.00308, 0.00329; %, 0.00330, 0.00330, 0.00330; %h-pol m1
%        0.00690, 0.00690, 0.00736, 0.00730, 0.00701; %, 0.00700, 0.00700, 0.00700; %v-pol m2
%        0.00600, 0.00600, 0.00656, 0.00660, 0.00660]'; %, 0.00660, 0.00660, 0.00660]'; %h-pol m2

%T_C=2.7;     %Background radiation [K]

%Emissivitet_FY_V=[0.9204; 0.9127; 0.9373; 0.9409; 0.9347]; % 0; 0; 0;];
%Emissivitet_FY_H=[0.7502; 0.7738; 0.8314; 0.8490; 0.8600]; % 0; 0; 0;];
%Emissivitet_MY_V=[0.9692; 0.9284; 0.8843; 0.8554; 0.7813]; % 0; 0; 0;];
%Emissivitet_MY_H=[0.8651; 0.8356; 0.7917; 0.7792; 0.7248]; % 0; 0; 0;];

%***********************
% Begin Calculation
%clear all
%close all

data_line=loadascii('/brain0/ltp/amsr/P1AME031118221MD_P01A0000000.00.ripTb.mask');

%data_line=loadascii('/brain2/amsr/orbits/2003.11.15/P1AME031115216MD_P01A0000000.00.ripTb')

%data_line=loadascii('C:\MATLAB6p5\work\AMSR\P1AME031027211MD.rip');
%load -ASCII C:\MATLAB6p5\work\Data\ice_orbit1.out
%load -ASCII ice_orbit1.out

%data_line=dlmread('C:\MATLAB6p5\work\Data\ice_orbit1.out',' '); 
%load data_line.mat

%data_line=ice_orbit1;

%Avv=[1.037, 1.032, 1.025, 1.032, 1.029];
%Ahv=[0.003, 0.003, 0.003, 0.004, 0.004];
%Aov=[-0.034, -0.029, -0.022, -0.028, -0.024];
%Ahh=[1.037, 1.031, 1.025, 1.034, 1.029];
%Avh=[0.003, 0.002, 0.003, 0.006, 0.004];
%Aoh=[-0.034, -0.029, -0.022, -0.028, -0.024];

len=size(data_line,1)
%T_A_V=[data_line(:,3),data_line(:,5), data_line(:,7), data_line(:,9), data_line(:,11)];
%T_A_H=[data_line(:,4),data_line(:,6), data_line(:,8), data_line(:,10), data_line(:,12)];

%for i=1:len
%    T_B_V(i,:)=(Avv.*T_A_V(i,:))+(Ahv.*T_A_H(i,:))+(2.7.*Aov);
%    T_B_H(i,:)=(Ahh.*T_A_H(i,:))+(Avh.*T_A_V(i,:))+(2.7.*Aoh);
%end

%data=[T_B_V(:,1), T_B_H(:,1), T_B_V(:,2), T_B_H(:,2), T_B_V(:,3), T_B_H(:,3), T_B_V(:,4), T_B_H(:,4), T_B_V(:,5), T_B_H(:,5)];
data=[data_line(:,3:12)];

%****************************************
% Start estimering

data_start=1;
data_slut=len;


T_A=[data(data_start:data_slut,:)];

[p, st, test]=inv_funktion2_apri_is(T_A(:,:));
%[p2, st2, fil_i2, fil_j2]=inv_funktion_apri_is_FY_MY(TA(:,:));
%[p1, st1]=inv_funktion_apri(T_A(:,:));

resultat=[data_line(data_start:data_slut,1:2), p(:,:), st(:,:), test(:,:)];
%resultat1=[data(data_start:data_slut,1:2), p1(:,:), st1(:,:)];

saveascii('/brain0/iomasa/dhb/Results/P1AME031118221MD_opti1.out',resultat)  


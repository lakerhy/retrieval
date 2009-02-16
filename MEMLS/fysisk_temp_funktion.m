%**************************************************************************
% FUNKTION DER BEREGNER DEN FYSISKE TEMPERATUR OG OVERFLADE TEMPERATUREN
% OVER EN ISOVERFLADE
%**************************************************************************
% Dorthe Hofman-Bang
% 12. november 2003
%**************************************************************************
% Metode: Iterativ metode baseret på Bootstrap algoritmen og emissiviteter
% Kilde: Algorithm Theoretical Basis Document (ATBD) for the AMSR_E Sea ice
%        algorithm; Donald J. Cavalieri and Josefino C. Comiso; Dec.2002
%**************************************************************************
% Input:  AMSR strålingstemperaturer
%         Format: T_B(6V, 6H, 10V, 10H, 18V, 18H, 23V, 23H, 37V, 37H)
% Output: Iskoncentration, fysisk temperatur og overfladetemperatur for en 
%         isoverflade
%         Format: [C_is, T_p, T_S]
%**************************************************************************
function [C_i3, C_i2, C_i1, T_p1, T_p2, T_S1, T_S2] = fysisk_temp_funktion(T_B)

%*******************************
% Bootstrap algorithm

TB6V=T_B(1,1);
%TB6H=T_B(1,2);
%TB10V=T_B(1,3);
%TB10H=T_B(1,4); 
TB18V=T_B(1,5);
%TB18H=T_B(1,6);  
%TB23V=T_B(1,7);  
%TB23H=T_B(1,8);  
TB37V=T_B(1,9);  
%TB37H=T_B(1,10);   

TW6V = 173.4639; %Beregnet fra data
TW37V = 208.9;   %218.8501 beregnet

TFY6V = 255.7488;    %beregnet fra data
TFY37V = 246.90; %255.7569

TMY6V =  260.6773;   %beregnet fra data
TMY37V = 185.40;

AF   = (TFY37V-TMY37V)/(TFY6V-TMY6V);                            
BF   = (TMY37V - AF*TMY6V);                                       
QF   = (TB37V-TW37V)/(TB6V-TW6V);                                
WF   = (TW37V - QF*TW6V);                                         
TI6VF = (BF-WF)/(QF-AF);                                          
TI37VF = AF*TI6VF+ BF;                                           
C_i1 = (TB6V - TW6V)/(TI6VF-TW6V);

%******************************************
% Calculate surface emissicity using a micing algorithm
% Emissivitet_FY_V=[0.9204; 0.9127; 0.9373; 0.9409; 0.9347; 0; 0; 0;];
% Emissivitet_FY_H=[0.7502; 0.7738; 0.8314; 0.8490; 0.8600; 0; 0; 0;];
% Emissivitet_MY_V=[0.9692; 0.9284; 0.8843; 0.8554; 0.7813; 0; 0; 0;];
% Emissivitet_MY_H=[0.8651; 0.8356; 0.7917; 0.7792; 0.7248; 0; 0; 0;];

epsilon6V_i=0.9549;  % Beregnet fra data
epsilon6V_0W=0.5937; % Beregnet fra data

epsilon1=epsilon6V_i*C_i1+epsilon6V_0W*(1-C_i1);

%******************************************
% Calculate the overall temperature og the surface within the pixel
T_p1=TB6V/epsilon1;

%******************************************
% Calculater eth emissicity og the surface at 18 and 37GHZ
epsilon_18V=TB18V/T_p1;
%epsilon_37H=TB37H/T_p1;
epsilon_37V=TB37V/T_p1;

%*****************************************
% Calculate ice concentraion, C_i using epsilon_18V, epsilon_37H and
% epsilon_37V and the Bootstrap Technique
e_FY18V=0.9373;
e_MY18V=0.8843;
e_W18V=0.6711;   % Beregnet fra data
e_FY37V=0.9347;
e_MY37V=0.7813;
e_W37V=0.7564;   % Beregnet fra data

AF   = (e_FY37V-e_MY37V)/(e_FY18V-e_MY18V);
BF   = (e_MY37V - AF*e_MY18V);                                       
QF   = (epsilon_37V-e_W37V)/(epsilon_18V-e_W18V);
WF   = (e_W37V - QF*e_W18V);                                         
TI18VF = (BF-WF)/(QF-AF);                                          
TI37VF = AF*TI18VF+ BF;                                           
C_i2 = (epsilon_18V - e_W18V)/(TI18VF-e_W18V);

T_S1=(T_p1-271*(1-C_i2))/C_i2;

%*****************************************************
epsilon2=epsilon6V_i*C_i2+epsilon6V_0W*(1-C_i2);

%******************************************
% Calculate the overall temperature og the surface within the pixel
T_p2=TB6V/epsilon2;

%******************************************
% Calculater eth emissicity og the surface at 18 and 37GHZ
epsilon2_18V=TB18V/T_p2;
%epsilon_37H=TB37H/T_p2;
epsilon2_37V=TB37V/T_p2;

%******************************************
% Calculate ice concentration using epsilon2_18V and epsilon2_37V

AF   = (e_FY37V-e_MY37V)/(e_FY18V-e_MY18V);
BF   = (e_MY37V - AF*e_MY18V);                                       
QF   = (epsilon2_37V-e_W37V)/(epsilon2_18V-e_W18V);
WF   = (e_W37V - QF*e_W18V);                                         
TI18VF = (BF-WF)/(QF-AF);                                          
TI37VF = AF*TI18VF+ BF;                                           
C_i3 = (epsilon2_18V - e_W18V)/(TI18VF-e_W18V);


%*****************************************
T_S2=(T_p2-271*(1-C_i3))/C_i3;

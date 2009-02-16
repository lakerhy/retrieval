%*************************************************************************
% FUNKTION TIL AT BEREGNE GEOFYSISKE PARAMETRE OVER IS OG ÅBENT VAND UD FRA 
% AMSR STRÅLINGSTEMPERATURER; A PRIORI INFORMATION ER INKLUDERET
%*************************************************************************
% Dorthe Hofman-Bang 
% 19. september 2003
%*************************************************************************
% Metode: Avanceret statistisk model for is og vand \m a priori information
% Kilde:  "Retrieval og Atmospheric Temperature and Composition From Remote 
%         Measurements of Thermal Radiation"; C.D.Rodgers; Reviews of
%         Geophysics and space Physics; Vol.14; NO.4; November 1976
%*************************************************************************
% Input:  Sæt SSMI strålingstemperaturer med formatet:
%         TA(6V,6H,10V,10H,18V, 18H, 23V, 23H, 37V, 37H)
% Output: fil_p:
%         Matrice med følgende parametre for hver SSMI måling
%               Vindhastighed, [m/s]
%               Integreret mængde af vanddamp i luftsøjlen, [mm]
%               Integreret mængde af flydende vand i luftsøjlen, [mm] 
%               T_ow Overfladetemperatur, [K]
%               T_is Overfladetemperatur, [K]
%               C_is iskoncentration i %
%               F_MY iskoncentration i %
%         fil_Std:
%           Matrice med standardafvigelserne for de 7 estimerede parametre
%           for hver af SSMI målingerne
%*************************************************************************
function [fil_p, fil_Std, fil_test] = inv_funktion2_apri_is(T_A)

teta=55*3.1415926535897/180;      % Indfaldsvinkel;  

lenT_A=length(T_A);

%*****************
%konstanter
%*****************
%Kovarians matricen af den normalfordelte fejl; Varians i diagonal
  %   [6V]    [6H]    [10V]   [10H]   [18V]   [18H]   [23V]   [23H]   [37V]   [37H]
S_e=[0.09,     0,      0,      0,      0       0,      0,      0,      0,      0;   
        0,     0.1089, 0,      0,      0       0,      0,      0,      0,      0;   
        0,     0,      0.2209, 0,      0       0,      0,      0,      0,      0;
        0,     0,      0,      0.2916, 0       0,      0,      0,      0,      0;
        0,     0,      0,      0,      0.2304  0,      0,      0,      0,      0;
        0,     0,      0,      0,      0,      0.2116, 0,      0,      0,      0;
        0,     0,      0,      0,      0,      0,      0.2025, 0,      0,      0;
        0,     0,      0,      0,      0,      0,      0,      0.1936, 0,      0;
        0,     0,      0,      0,      0,      0,      0,      0,      0.2025, 0;  
        0,     0,      0,      0,      0,      0,      0,      0,      0,      0.16;];

%*****************
% A priori information
%Kovarians matrice (middelværdier er givet  inde i løkken)

S_p=[ 12.3024,    3.2072,    0.1398,      6.0322,          0,     -0.6525,  -0.9347;
       3.2072,   11.0481,    0.2495,     11.9348,          0,     -0.3085,  -0.5362;
       0.1398,    0.2495,    0.0204,      0.2041,          0,     -0.0063,  -0.0152;
       6.0322,   11.9348,    0.2041,     23.9468,          0,     -0.7254,  -1.0207;
            0,          0,         0,           0,    23.9468,          0,         0;
      -0.6525,   -0.3085,   -0.0063,     -0.7254,          0,      0.0114,   0.1203;
      -0.9347,   -0.5362,   -0.0152,     -1.0207,          0,      0.1203,   0.0332];
S_p_inv=S_p^(-1);
S_e_inv=S_e^(-1);
%***************************
%Define matix

fil_p=zeros(lenT_A,7);
fil_Std=zeros(lenT_A,7);
fil_test=zeros(lenT_A,1);
%****************************
% Start beregninger
%****************************    
for nn=1:lenT_A

    if ( mod(nn,1000) == 0 ) %Udskriv hver 100 nr. til screen
        nn  
    end

    %*******************************
    % Beregn startgæt
    %*******************************
    % Beregn startgæt på iskoncentrationer
    [C_FY1, C_MY1, C_is1]=is_C_funktion_AMSR(T_A(nn,:));  %Beregn iskoncentration for pågældende måling (IKKE %)
       
    if C_is1<=0.001
        C_is1=0.001;
        F_MY=0.001;
    end
        
    if C_MY1<=0.001
        F_MY1=0.001;
    end
    
    if C_MY1>=C_is1
        if C_is1<0.1
            F_MY1=0.001;
        else
            F_MY1=C_is1;
        end
    else
        if C_is1<0.1
            F_MY1=0.001;
        else
            F_MY1=C_MY1/C_is1;
        end
    end
    
    %Beregn startgæt på T_is
    [C_i3, C_i2, C_i1, T_p1, T_p2, T_S1, T_S2] = fysisk_temp_funktion(T_A(nn,:));
    
    if C_i2>0.5
        T_is=T_S2;    
    else
        T_is=273;
    end
    
    % Startgæt, 
    p=[        4.75,          5.25,      5.0;    %W, vindhastighed   5; 10;  8
               3.40,          3.80,     3.60;    %V, vanddamp        2; 8;   4
             0.0775,        0.0825,     0.08;    %L, Liquid water    0; 0.2; 0     
                274,           275,    274.5;    %T_s, overflade temperatur 
         (T_is-0.5),    (T_is+0.5),     T_is;    %T_is, is temperatur    
       (C_is1-0.01),  (C_is1+0.01),    C_is1;    %C_FY, Koncentration af FY is
       (F_MY1-0.01),  (F_MY1+0.01),    F_MY1;];  %C_MY, Koncentration af MY is
    
    % A priori middelværdier (for isen NASA Team algoritmen)   
    p_0=[   4.9533;    3.6164;    0.0808;  274.5;  T_is;   C_is1;    F_MY1;];
 
    %**************************************************
    %Beregn T_A0 for gæt på p0 v.h.a forward algoritme; 
    %**************************************************

    [T_A0] = FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)); %W,V,L,T_s,C_is,F_MY; kald til forward funktion
   
    %*********************************
    %Beregn fejlen mellem T_A og T_A0
    %*********************************
   
    Delta_T=T_A(nn,:)'-T_A0;              % Find fejlen
   
    %*********************************
    %Start beregning af nyt p estimat 
    %*********************************
    for ite=1:5
    	
       	%Bergen M matrice
       	M(:,1)=(FW_funktion2_is(p(1,1),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3))-FW_funktion2_is(p(1,2),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(1,1)-p(1,2));  %dT/dW
       	M(:,2)=(FW_funktion2_is(p(1,3),p(2,1),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3))-FW_funktion2_is(p(1,3),p(2,2),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(2,1)-p(2,2));  %dT/dV
       	M(:,3)=(FW_funktion2_is(p(1,3),p(2,3),p(3,1),p(4,3),p(5,3),p(6,3),p(7,3))-FW_funktion2_is(p(1,3),p(2,3),p(3,2),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(3,1)-p(3,2));  %dT/dL
       	M(:,4)=(FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,1),p(5,3),p(6,3),p(7,3))-FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,2),p(5,3),p(6,3),p(7,3)))/(p(4,1)-p(4,2));  %dT/dT_s
       	M(:,5)=(FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,1),p(6,3),p(7,3))-FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,2),p(6,3),p(7,3)))/(p(5,1)-p(5,2));  %dT/dC_is
       	M(:,6)=(FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,1),p(7,3))-FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,2),p(7,3)))/(p(6,1)-p(6,2));  %dT/dC_FY
       	M(:,7)=(FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,1))-FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,2)))/(p(7,1)-p(7,2));  %dT/dC_MY
       
        Delta_p=p_0-p(:,3); %Beregn fejlen mellem p_0 og p
       
        %***************
        %Beregn delta p
        %***************
                
        S=((S_p_inv)+M'*(S_e_inv)*M)^-1;                           % Kovariansmartrice
        p_est=S*(((S_p_inv)*Delta_p)+(M'*(S_e_inv)*Delta_T));      % LS estimat af vektor med goefysiske parametre
                
        %********************************************************
        % Addder delta p til p0 med mindre min  værdien er naaet
        %********************************************************
        for par=1:3     %Kør gennem de 3 parametre
            if p(par,2)+p_est(par)<=0           %Alle parametre er mindre en nul 
                p(par,:)=[0.001, 0.003, 0.002]; %Min, Max, Default
            elseif p(par,3)+p_est(par)<=0       %Default og min er mindre end nul
                p(par,1)=0.001; %Min
                p(par,3)=0.002; %Default
                                     %Max
                if p(par,2)+p_est(par)>p(par,3) %Tjek at den nye maxværdi vil blive større end defaultværdien 
               	    p(par,2)=p(par,2)+p_est(par);
                else
                    p(par,2)=0.003;
                end                  
            elseif p(par,1)+p_est(par)<=0       %Minværdien er mindre end nul
                p(par,1)=0.001;  %Min hæves
                                    %Max og Default
                if p(par,3)+p_est(par)>p(par,1) %Tjek at minværdien er mindre end den nye defaultværdien     
               	    p(par,3)=p(par,3)+p_est(par);   %Default
                    p(par,2)=p(par,2)+p_est(par);   %Max
                elseif p(par,3)+p_est(par)<p(par,1) & p(par,1)<p(par,2)+p_est(par)%Tjek at ny min er større end den ny default og mindre end ny max     
                    p(par,3)=0.002; %Default
                    if p(par,3)<(p(par,2)+p_est(par))   %Tjek at den ny defaultværdi er mindre end max
                  	    p(par,2)=p(par,2)+p_est(par);%Max
                    else
                        p(par,2)=0.003;
                    end
                else                                %Ny Minværdi er større end ny Maxværdi
                    p(par,3)=0.002;
                    p(par,2)=0.003;
                end
            else %Min større end nul
                p(par,:)=p(par,:)+p_est(par);
            end   
        end
        for par=4:7
            if par==4  %T_S, overflade temp
                %if p(par,3)+p_est(par)>=271   %T_S skal være større end -2deg C    
                    p(par,:)=p(par,:)+p_est(par);       %Normal opdatering
                    %else
                    %p(par,:)=[269.9, 271.1, 271.0];     %Mingrænse er nået
                    %end
            elseif par==5   %T_is, Is temperatur
                %if p(par,3)+p_est(par)<273   %T_is skal være mindre end 0 [deg C]
                    p(par,:)=p(par,:)+p_est(par);       %Normal opdatering
                    %else
                   % p(par,:)=[272.9, 273.1, 273];       %Maksgrænse er nået
                   %end
            elseif par==6    %C_is, iskoncentration
                 p(par,:)=p(par,:)+p_est(par);           %Normal opdatering
            elseif par==7    % F_MY, FY ice fraction
                 if (p(par,3)+p_est(par))>1
                    p(par,:)=[0.9991, 0.9993, 0.9992];  %Maks grænse er nået
                elseif (p(par,3)+p_est(par))<0
                    p(par,:)=[0.0001, 0.0003, 0.0002];  %Min grænse er nået
                else
                    p(par,:)=p(par,:)+p_est(par);   %Normal opdatering
                end
            end
        end
        %*****************************************************
        %Beregn T_A0 for nye p værdier v.h.a forward algoritme 
        %*****************************************************
        [T_A0] = FW_funktion2_is(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3));     % kald til forward funktion
        
        %**********************************************
        %Beregn fejlen mellem T_A og T_A0
        %**********************************************
        Delta_T=T_A(nn,:)'-T_A0;              % Find fejlen

        ite_Std(ite,:)=[sqrt(S(1,1)), sqrt(S(2,2)), sqrt(S(3,3)), sqrt(S(4,4)), sqrt(S(5,5)), sqrt(S(6,6)), sqrt(S(7,7))];%Gem variansen for hvert datasæt [W,V,L,T_s,C_FY,C_MY]
        %ite_Delta_T(ite,:)=Delta_T';    %Test variabel
        ite_p(ite,:)=p(:,3)';           %opsaml p for hvert enkelt intration
        %ite_pest(ite,:)=p_est';        %Test variabel
        %ite_TA0(ite,:)=T_A0;           %Test variabel
        
       test(ite)=sqrt(sum(Delta_T.^2));
      %test(ite)=sum(abs(Delta_T));
    end
    
   % FMY(nn,:)=F_MY1;            %Gem variabel fra hver T_B med F_MY
   % CIS(nn,:)=C_is1;            %Gem Variabel fra  hver T_B med hver F_MY
   
    [i,j]=min(test);            %Gem det estimat med den mindste fejl
    fil_test(nn,:)=i;              %returner testvariabel
   % fil_j(nn,:)=j;              %returner indeks på testvariabel
    
   % test_all(nn,:)=test;
    
    fil_p(nn,:)=ite_p(j,:);     %opsaml p for hvert enkelt datasæt
    fil_Std(nn,:)=ite_Std(j,:); %opsaml std for hvert enkelt datasæt
  
end

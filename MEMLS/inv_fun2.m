% Input:  The measured radiometer tempetature vector TA:
%         The elements of the vector is made up accordingly
%         TA(6V, 6H, 10V, 10H, 18V, 18H, 23V, 23H, 37V, 37H)
%          where 6V is the vertical polarisation of the 6GHz channel. 
% Output: Geophysical parameters estimated from input temperature
%               T_s, surface Temperature [K]
%               W,water content of air [mm]
%              di_snow, Thickness of snow [cm]
%               roi_snow, Snow density [kg/m2]
%               pci_snow, Snow correlation length [mm] 
%               di_ice, Thickness of ice [cm]
%               sal [per mil]
% 
%************************************************************************* 
function [p_estimation,TraceS] = inv_fun2(T_A) 

%Covariance matrix of the AMSR measurements
%   [6V]    [6H]    [10V]   [10H]   [18V]   [18H]   [23V]   [23H]   [37V]   [37H]
  S_e=[0.09,     0,      0,      0,      0       0,      0,      0,      0,      0;  
       0,     0.1089, 0,      0,      0      0,      0,      0,      0,      0;  
       0,     0,      0.2209, 0,      0       0,      0,      0,      0,      0;
       0,     0,      0,      0.2916, 0       0,      0,      0,      0,      0;
       0,     0,      0,      0,      0.2304  0,      0,      0,      0,      0;
       0,     0,      0,      0,      0,      0.2116, 0,      0,      0,      0;
       0,     0,      0,      0,      0,      0,      0.2025, 0,      0,      0;
       0,     0,      0,      0,      0,      0,      0,      0.1936, 0,      0;
       0,     0,      0,      0,      0,      0,      0,      0,      0.2025, 0; 
       0,     0,      0,      0,      0,      0,      0,      0,      0,      0.16;]; 
  %Exclude channel
  % To exclude a channel uncomment the appropriate S_e index 
  %S_e(1,1)=10000;    % 6V
  %S_e(2,2)=10000;    % 6H
  %S_e(3,3)=10000;     %10V
  %S_e(4,4)=10000;
  %S_e(5,5)=10000;     %18V
  %S_e(6,6)=10000;
  %S_e(7,7)=10000;     %23V
  %S_e(8,8)=10000;
  %S_e(9,9)=10000;     %37V
  %S_e(10,10)=10000; 
  %*****************
  % A priori information
  %Covariance matrix for the a priori information 
  % S_p=[ 12.3024,    3.2072,    0.1398,      6.0322,          0,     -0.6525,  -0.9347;
  %        3.2072,   11.0481,    0.2495,     11.9348,          0,     -0.3085,  -0.5362;
  %        0.1398,    0.2495,    0.0204,      0.2041,          0,     -0.0063,  -0.0152;
  %        6.0322,   11.9348,    0.2041,     23.9468,          0,     -0.7254,  -1.0207;
  %             0,          0,         0,           0,    23.9468,          0,         0;
  %       -0.6525,   -0.3085,   -0.0063,     -0.7254,          0,      0.0114,   0.1203;
  % %       -0.9347,   -0.5362,   -0.0152,     -1.0207,          0,      0.1203,   0.0332];

  load dev_corr;

  S_p = covariance(dev_corr);

  S_p_inv=S_p^(-1);
  S_e_inv=S_e^(-1); 
  %*************************************************************************
  %Find Iteration starting point
  %**************************************************************************
  %First guess

  p=[245,      246,      245.5;                  % T_s
       1.800,    1.810,    1.805;                  % W 
       22,       23,       22.5;                   % Di_snow              
       346,      347,      346.5;                  % Roi_snow
       0.13,     0.135,    0.1325;                 % Pci_snow
       477,      478,      477.5;                  % Di_ice
       3.3,      3.31,     3.305;]                 % Sal
  
  % A priori matrix
  p
  p_0=[0.2456533956128566e3,0.1892598074880626e1,0.2241547871754426e2, ...
       0.3462067506803534e3,0.1294358595730483, 0.4775764424657535e3,0.3348532459662817e1]';

  %*****************************************************************************
  % Calculate T_A_0    Temperature for guessed value
  % and Delta_T           The error of the first guess
  %**************************************************************************
  T_A_0 = fw_fun2(p_0(1,1),p(2,1),p_0(3,1),p_0(4,1),p_0(5,1),p_0(6,1),p_0(7,1));
  Delta_T=T_A-T_A_0; 
  
  if max(abs(Delta_T)) < 0.5
    iterate=0;
    p_est=p_0;             %Return with p
  else
    iterate=1;
  end
  k=0;
  a=0;

  %*********************************************************************************
  %       ITERATION
  %*********************************************************************************
  while iterate
    k=k+1;
    a=a+1;
    %**********************************************************************
    %    Calculate M    Mixing matrix for p_iteration
    %********************************************************************************
    %steps in the calculation of the dT/dW
    
    % DT/DT_s
    M(:,1)=(fw_fun2(p(1,1),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3))-fw_fun2(p(1,2),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(1,1)-p(1,2));
    % DT/DW        
    M(:,2)=(fw_fun2(p(1,1),p(2,1),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3))-fw_fun2(p(1,2),p(2,2),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(2,1)-p(2,2));
    % DT/DDi_snow
    M(:,3)=(fw_fun2(p(1,1),p(2,3),p(3,1),p(4,3),p(5,3),p(6,3),p(7,3))-fw_fun2(p(1,2),p(2,3),p(3,2),p(4,3),p(5,3),p(6,3),p(7,3)))/(p(3,1)-p(3,2));
    % DT/DRoi_snow    
    M(:,4)=(fw_fun2(p(1,1),p(2,3),p(3,3),p(4,1),p(5,3),p(6,3),p(7,3))-fw_fun2(p(1,2),p(2,3),p(3,3),p(4,2),p(5,3),p(6,3),p(7,3)))/(p(4,1)-p(4,2));
    % DT/DPci_snow    
    M(:,5)=(fw_fun2(p(1,1),p(2,3),p(3,3),p(4,3),p(5,1),p(6,3),p(7,3))-fw_fun2(p(1,2),p(2,3),p(3,3),p(4,3),p(5,2),p(6,3),p(7,3)))/(p(5,1)-p(5,2));
    % DT/DDi_ice    
    M(:,6)=(fw_fun2(p(1,1),p(2,3),p(3,3),p(4,3),p(5,3),p(6,1),p(7,3))-fw_fun2(p(1,2),p(2,3),p(3,3),p(4,3),p(5,3),p(6,2),p(7,3)))/(p(6,1)-p(6,2));
    % DT/Dsal    
    M(:,7)=(fw_fun2(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,1))-fw_fun2(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,2)))/(p(7,1)-p(7,2));
    %**********************************************************************
    %    Calculate delta_p    Calculate what to add to p-iteration to get
    %    to p-true for TA with the current mixing matrix
    %********************************************************************************
    delta_p=p_0-p(:,3); 
    %********************************************************************************
    %    Calculate next iteration p   
    %******************************************************************************** 
    S=((S_p_inv)+M'*(S_e_inv)*M)^-1;
    p_est = S*(((S_p_inv)*delta_p)+(M'*(S_e_inv)*Delta_T'));

    TraceS=trace(abs(S));
    
    % validality conditions for updating
    for par=1:7
      
      % T_s [229.95,271.15]
      if par==1
        if p(par,3)+p_est(par) >= 271.15
          p(par,:)=[271.05,271.15,271.1];
        elseif p(par,3)+p_est(par) <= 229.95
          p(par,:)=[229.95,230.05,229.90];
        else p(par,:)=p(par,:)+p_est(par);
        end
      end
      
      if par==2|par==7|par==5
        if p(par,2)+p_est(par)<=0
          p(par,:)=[0.001,0.003,0.002];
        elseif p(par,3)+p_est(par)<=0
          p(par,1)=0.001;
          p(par,3)=0.002;
          if p(par,2)+p_est(par)>p(par,3)
            p(par,2)=p(par,2)+p_est(par);
          else 
            p(par,2) = 0.003;
          end
        elseif p(par,1)+p_est(par)<=0
          p(par,1)=0.001;
          if p(par,3)+p_est(par)>p(par,1)
            p(par,3)=par(par,3)+p_est(par);
            p(par,2)=par(par,2)+p_est(par);
          elseif p(par,3)+p_est(par)<p(par,1) & p(par,1)<p(par,2)+ ...
                p_est(par)
            p(par,3)=0.002;
            if p(par,3)<(p(par,2)+p_est(par))
              p(par,2)=p(par,2)+p_est(par);
            else p(par,2) = 0.002;
            end
          else
            p(par,3)=0.002;
            p(par,2)=0.003;
          end
        else
          p(par,:) = p(par,:)+p_est(par);
        end
      end
      
      if par==3|par==4|par==6
        if p(par,2)+p_est(par)<=0
          p(par,:)=[1,3,2];
        elseif p(par,3)+p_est(par)<=0
          p(par,1)=1;
          p(par,3)=2;
          if p(par,2)+p_est(par)>p(par,3)
            p(par,2)=p(par,2)+p_est(par);
          else 
            p(par,2) = 3;
          end
        elseif p(par,1)+p_est(par)<=0
          p(par,1)=0.1;
          if p(par,3)+p_est(par)>p(par,1)
            p(par,3)=p(par,3)+p_est(par);
            p(par,2)=p(par,2)+p_est(par);
          elseif p(par,3)+p_est(par)<p(par,1) & p(par,1)<p(par,2)+ ...
                p_est(par)
            p(par,3)=2;
            if p(par,3)<(p(par,2)+p_est(par))
              p(par,2)=p(par,2)+p_est(par);
            else p(par,2) = 2;
            end
          else
            p(par,3)=2;
            p(par,2)=3;
          end
        else
          p(par,:) = p(par,:)+p_est(par);
        end
      end

%       % pci_snow [-1,1]
%       if par==5
%         if p(par,3)+p_est(par)< -1
%           p(par,:)=[-0.999,-0.997,-0.998];
%         elseif p(par,3)+p_est(par)> 1
%           p(par,:)=[0.997,0.999,0.998];
%         else 
%           p(par,:) = p(par,:)+p_est(par);
%         end
%       end
      
    end
    
    T_A_estimated=fw_fun2(p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3));
    Delta_T=(T_A-T_A_estimated)
    
    iteration_std(k,:)=[sqrt(S(1,1)),sqrt(S(2,2)),sqrt(S(3,3)),sqrt(S(4, ...
                                                      4)),sqrt(S(5,5)),sqrt(S(6,6)),sqrt(S(7,7))]
    T_std(k)=sqrt(sum(Delta_T.^2))
    
    %Exit conditions could be set either to depend on Delta_T or be a fixed
    %number of times.
    if max(abs(Delta_T)) < 0.1
      iterate=0;
      p_estimation = [p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)];
    else
      iterate=1;
    end
    if a==3
      iterate=0;
            p_estimation = [p(1,3),p(2,3),p(3,3),p(4,3),p(5,3),p(6,3),p(7,3)];
    end
  end 
  


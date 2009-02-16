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
function [p,TraceS] = inv_fun(T_A,type) 

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

  
  if type ==1 %% MY
    load dev_corr_MY;

    S_p = covariance(dev_corr_MY);
    p=[ 245.5;                  % T_s
        1.805;                  % W 
        22.5;                   % Di_snow              
        346.5;                  % Roi_snow
        0.1325;                 % Pci_snow
        477.5;                  % Di_ice
        3.305;]                 % Sal
    
    p0=[0.2456533956128566e3,0.1892598074880626e1,0.2241547871754426e2, ...
         0.3462067506803534e3,0.1294358595730483, 0.4775764424657535e3, ...
         0.3348532459662817e1]';  %%%% For MY start mean value

  end
  
  if type ==2 %% FY
    load dev_corr_FY;

    S_p = covariance(dev_corr_FY);
    p=[ 245.5;                  % T_s
        1.805;                  % W 
        22.5;                   % Di_snow              
        346.5;                  % Roi_snow
        0.2125;                 % Pci_snow
        140;                  % Di_ice
        6.705;]                 % Sal
    p0=[0.2456533956128566E+03,0.1892598074880626E+01, 0.2153510547453996E+02, ...
         0.3447630980460535E+03, 0.2125703066295413E+00, 0.1412827332719277E+03, ...
         0.6719991046581256E+01]'; %% For FY start mean value

  end

  S_p_inv=S_p^(-1);
  S_e_inv=S_e^(-1); 
  %*************************************************************************
  %Find Iteration starting point
  %**************************************************************************
   %*****************************************************************************
  % Calculate T_A_0    Temperature for guessed value
  % and Delta_T           The error of the first guess
  %**************************************************************************
  T_A_0 = fw_fun3(p,type);
  Delta_T=T_A-T_A_0; 
  if max(abs(Delta_T)) < 0.5
    iterate=0;
    else
    iterate=1;
                 %Start iteration with p
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
    dT        = [.1/k;0;0;0;0;0;0];
    dW        = [0;0.05/k;0;0;0;0;0];
    ddi_snow  = [0;0;.1/k;0;0;0;0];
    droi_snow = [0;0;0;.1/k;0;0;0];
    dpci_snow = [0;0;0;0;0.001/k;0;0];
    ddi_ice   = [0;0;0;0;0;.1/k;0];
    dsal      = [0;0;0;0;0;0;0.01/k];
    flag_k=k

    M(:,1)=(fw_fun3(p+dT/2,type)-fw_fun3(p-dT/2,type))/dT(1,1);         
    M(:,2)=(fw_fun3(p+dW/2,type)-fw_fun3(p-dW/2,type))/dW(2,1);         
    M(:,3)=(fw_fun3(p+ddi_snow/2,type)-fw_fun3(p-ddi_snow/2,type))/ddi_snow(3,1);          
    M(:,4)=(fw_fun3(p+droi_snow/2,type)-fw_fun3(p-droi_snow/2,type))/droi_snow(4,1);     
    M(:,5)=(fw_fun3(p+dpci_snow/2,type)-fw_fun3(p-dpci_snow/2,type))/dpci_snow(5,1);  
    M(:,6)=(fw_fun3(p+ddi_ice/2,type)-fw_fun3(p-ddi_ice/2,type))/ddi_ice(6,1); 
    M(:,7)=(fw_fun3(p+dsal/2,type)-fw_fun3(p-dsal/2,type))/dsal(7,1);    
    %**********************************************************************
    %    Calculate delta_p    Calculate what to add to p-iteration to get
    %    to p-true for TA with the current mixing matrix
    %********************************************************************************
    delta_p=p0-p; 
    %********************************************************************************
    %    Calculate next iteration p   
    %******************************************************************************** 
    S=((S_p_inv)+M'*(S_e_inv)*M)^-1;
    p=p+S*(((S_p_inv)*delta_p)+(M'*(S_e_inv)*Delta_T'));

    TraceS=trace(abs(S));
    
    %     if p(1)< 240
    %       p(1) = 240+1/k;
    %     end
    if p(1) > 271.15 
      p(1)=271.15;
    elseif p(1) < 229.95
      p(1)=229.95;
    end
    
    if p(2)< 0
      p(2) = 0.5/k;
    end
    if p(3)<0
      p(3) = 1/k;
    end
    
    if p(4)< 0
      p(4) = 1/k;
    end
    
    if p(5)<-1
      p(5) = -1+0.005/k;
    elseif p(5)>1
      p(5) = 1-0.005/k;
    end
    
    if p(6)< 0
      p(6) = 1/k;
    end
    
    if p(7)< 0
      p(7) = 0.5/k;
    end

    T_A_estimated=fw_fun3(p,type);
    
    Delta_T=(T_A-T_A_estimated)
    
    iteration_std(k,:)=[sqrt(S(1,1)),sqrt(S(2,2)),sqrt(S(3,3)),sqrt(S(4, ...
                                                      4)),sqrt(S(5,5)),sqrt(S(6,6)),sqrt(S(7,7))]
    T_std(k)=sqrt(sum(Delta_T.^2))
    
    
    %Exit conditions could be set either to depend on Delta_T or be a fixed
    %number of times.
    if max(abs(Delta_T)) < 1
      iterate=0;
    else
      iterate=1;
    end
    if a==1
      iterate=0;
    end
  end 
  


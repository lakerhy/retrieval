
% Input:  The measured radiometer tempetature vector TA:
%         The elements of the vector is made up accordingly   TA(6V, 6H, 10V, 10H, 18V, 18H, 23V, 23H, 37V, 37H)
%          where 6V is the vertical polarisation of the 6GHz channel. 
% Output: Geophysical parameters estimated from input temperature
%               T_s, surface Temperature [K]
%              di_snow, Thickness of snow [cm]
%               pci_snow, Snow correlation length [mm] 
%               roi_snow, Snow density [kg/m3]
%               CFY
%************************************************************************* 
function [p_estimation,std_S,std_Sp] = inversion(T_A) 

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
  S_p = zeros(5,5);
  S_p(1,1) = 100;
  S_p(2,2) = 100;
  S_p(3,3) = 0.005;
  S_p(4,4) = 1600;
  S_p(5,5) = 0.25;
  diagS_p = diag(S_p);
  std_Sp= sqrt(diagS_p);
  S_p_inv=S_p^(-1);
  S_e_inv=S_e^(-1); 
  %*************************************************************************
  %Find Iteration starting point
  %**************************************************************************
  %First guess
  %   p=[249,      250,      249.5;                  % T_s
  %       17.9,    18,    17.95;       
  %        0.09,       0.092,       0.091;  
  %        319.8,      320,      319.9;            
  %        0.498,     0.5,    0.499;]
  % A priori matrix
  p0=[235,18,0.1,320,0.1]';
  p=[231,30,0.11,325,0]';
  %*****************************************************************************
  % Calculate T_A_0    Temperature for guessed value
  % and Delta_T           The error of the first guess
  %**************************************************************************
  T_A_0 = fw(p);

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
    
    dT        = [0.1/k;0;0;0;0];
    ddi_snow  = [0;.05/k;0;0;0];
    dpci_snow = [0;0;0.0005/k;0;0];
    droi_snow = [0;0;0;0.5/k;0];
    dCFY   = [0;0;0;0;.0005/k];
    

    M(:,1)=(fw(p+dT/2)-fw(p-dT/2))/dT(1,1);
    
    M(:,2)=(fw(p+ddi_snow/2)- fw(p-ddi_snow/2))/ddi_snow(2,1);
    M(:,3)=(fw(p+dpci_snow/2)-fw(p-dpci_snow/2))/dpci_snow(3,1) ;
    M(:,4)=(fw(p+droi_snow/2)-fw(p-droi_snow/2))/droi_snow(4,1)     ;
    M(:,5)=(fw(p+dCFY/2)-fw(p-dCFY/2))/dCFY(5,1)  ;
    
    %**********************************************************************
    %    Calculate delta_p    Calculate what to add to p-iteration to get
    %    to p-true for TA with the current mixing matrix
    %********************************************************************************
    delta_p=p0-p; 
    %********************************************************************************
    %    Calculate next iteration p   
    %******************************************************************************** 
    S=((S_p_inv)+M'*(S_e_inv)*M)^-1;
    
    p = p+S*(((S_p_inv)*delta_p)+(M'*(S_e_inv)*Delta_T'))
     diagS=diag(S);
     std_S = sqrt(diagS); 
    TraceS=trace(abs(S));
    
    % validality conditions for updating
    
    
    
    if p(1) > 271.15 
      p(1)=271.15;
    elseif p(1) < 229.95
      p(1)=229.95;
    end
    
    if p(2)< 0
      p(2) = 0.1/k;
    end
    if p(3)<0
      p(3) = 0.01/k;
    end
    
    if p(4)< 0
      p(4) = 1/k;
    end
    
    if p(5)<0
      p(5) = 0.001/k;
    elseif p(5)>1
      p(5) = 1-0.001/k;
    end
    
    
    T_A_estimated=fw(p);
    Delta_T=(T_A-T_A_estimated);
    
    iteration_std(k,:)=[sqrt(S(1,1)),sqrt(S(2,2)),sqrt(S(3,3)),sqrt(S(4,4)),sqrt(S(5,5))];
    T_std(k)=sqrt(sum(Delta_T.^2));
    
    %Exit conditions could be set either to depend on Delta_T or be a fixed
    %number of times.
    if max(abs(Delta_T)) < 1
      iterate=0;
      p_estimation = p;
    else
      iterate=1;
    end
    if a==6
      iterate=0;
      p_estimation = p;
        end
  end 
  
  
  


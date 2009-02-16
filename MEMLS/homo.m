% ************************* Usage *****************************
% output the memls layer profile file 
% ************************* Input *****************************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               [K]        [kg/m^-3]    [cm]     [mm] [per mil] ice/snow[1/0]  
%%         num,   T     W   ,  roi   ,    di ,     pci   sal     type
%%  par = [1 , Tsea ,    *, roi_ice,   di_ice,   *,      7.5,     1;
%          2 , *    , W_slush,roi_slush,di_slush,pci_slush,     0  1
%%         3 , Ts_snow , 0, roi_snow, di_snow, pci_snow, 0,       0;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snow_layernums: numbers of snow layers
% *************************Output*****************************
% syntes.disc: output_file with layer- number, temp [K], vol.water content, density [kg/m3], thickness [cm], cor.length [mm],  salinity [per mil], ice/snow[1/0]
function[par_out] = homo (par_in,snow_layernums,ice_type)

%%%%%%%%%%%%%%CONSTANTS%%%%%%%%%%%%%%%%%%%%%%%
  Tsea=271.35;      %[K] 
  ks=0.31;
  k0=2.034;         %[W/m/K]
  beta=0.117        %[W/m/per mil] 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%  
  Ts_snow = par_in(3,2); 
  roi_snow = par_in(3,4); 
  di_ice = par_in(1,5);  
  di_snow = par_in(3,5);
  di_slush=par_in(2,5);
  pci_snow = par_in(3,6); 
  sal_ice = par_in(1,7); 
  W_slush = par_in(2,3);
  roi_slush = par_in(2,4);
  pci_slush = par_in(2,6);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % Range of validity of surface tempurature of snow surface
  if Ts_snow > 271.15
    Ts_snow = 271.35;
  elseif Ts_snow < 229.95
    Ts_snow = 229.95;
  end

  T(1) = Tsea;
  for i = 1:80
    %    ki = k0 + (beta*sal_ice)/(T(1)+273.15); %conductivity of ice 
    ki = k0 + (beta*sal_ice)/T(1); %conductivity of ice 
    slope = ki*di_ice/ks/(di_snow+di_slush); %Ratio between the temperature gradients in ice and snow
    T0(1) = (Ts_snow + slope*Tsea)/(1 + slope); %Temperature of the interface of
                                                %snow and ice
    slope2 = (T0(1)-Tsea)/di_ice; % temperature gradient of ice layer
    if abs(slope2) > 0.1 %4.15f
      slope2 = slope2/abs(slope2)*0.1;
      T0(1) = Tsea+slope2*di_ice;
    end
    slope3 = (Ts_snow-T0(1))/(di_snow+di_slush);  % temperture slope of snow layer
    T(1) = Tsea+slope2*di_ice/2;
    T(2) = T0(1)+slope3*(di_snow+di_slush)/2;
    
    T(1) = (Tsea+T(1))/2;  % mean value of ice temperature
                           %T0(2)=T0(1)
    if (i>1 && abs(T0(1)-T0(2))<0.01) 
      break
    end
    T0(2)=T0(1);
    
  end 
  
  T=T0(1);  % interface temperature
  T_ice=T(1); % ice temperature
  T_slush = (di_snow*T+Ts_snow*di_slush)/(di_snow+di_slush);
  delta_T1 = (T_slush - Ts_snow)/snow_layernums; 
  delta_T2 = Tsea - T;  
  

  Ti(1) = Tsea - delta_T2/2; % average ice temperature
  
  Ti(2)= (T_slush+T)/2; % average slush temperature

  for i = 1:(snow_layernums+1)
    Ti(i+2) = T_slush - (i-1)*delta_T1;
  end
  for i = 1:(snow_layernums)
    Ti(i+2) = (Ti(i+2) + Ti(i+3))/2;
  end
 
  T_ice = Ti(1);
 
  
  T0 = abs(T_ice - 273.15);
  if T0<1
    W_ice = 1e-3*sal_ice*49.717;
  else 
    W_ice = 1e-3*sal_ice*(-49.185/T_ice+0.532);
  end
  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  layer_num = [1:snow_layernums+2];
  W = [W_ice,W_slush,zeros(1,snow_layernums)];  %Water conetent of ice

  if ice_type==1 %MY
    if di_ice >= 25
      pci_MY = 0.35;
    else 
      pci_MY = (0.35*0.25+0.25*(di_ice/100-0.25))*100/di_ice;
    end
    roi_MY = 918;

    roi = [roi_MY,roi_slush,roi_snow*ones(1,snow_layernums)];
    di = [di_ice,di_slush,di_snow/snow_layernums*ones(1,snow_layernums)];
    pci = [pci_MY,pci_slush,pci_snow*ones(1,snow_layernums)];
    sal = [sal_ice,0,zeros(1,snow_layernums)];
    si = [1,0,zeros(1,snow_layernums)];

    par_out = [layer_num;Ti(1:snow_layernums+2);W;roi;di;pci;sal;si];
    
    %     fid = fopen(output_file,'wt');
    %     output = [layer_num;Ti(1:snow_layernums+1);W;roi;di;pci;sal;si];
    %     fprintf(fid,'%3d %3.12f %3.12f %4.12f %4.15f %4.15f %3.15f %2d\n', output);
    %     fclose(fid);
    
    %     fid = fopen('input_MY', 'wt');
    %     MY_input = [Tsea,Ts_snow;roi_MY,roi_snow;di_ice,di_snow;pci_MY,pci_snow;sal_ice,0];
    %     fprintf(fid,'%3.12f %3.12f %4.12f %4.15f %4.15f\n',MY_input);
    %     fclose(fid);

  end

  if ice_type==2 %FY
    
    pci_FY = 0.8*(3*W_ice/(4*3.14))^0.333;
    roi_FY = 920;
    
    roi = [roi_FY,roi_slush,roi_snow*ones(1,snow_layernums)];
    di = [di_ice,di_slush,di_snow/snow_layernums*ones(1,snow_layernums)];
    pci = [pci_FY,pci_slush,pci_snow*ones(1,snow_layernums)];
    sal = [sal_ice,0,zeros(1,snow_layernums)];
    si = [1,0,zeros(1,snow_layernums)];

    
    par_out = [layer_num;Ti(1:snow_layernums+2);W;roi;di;pci;sal;si];
    
  end

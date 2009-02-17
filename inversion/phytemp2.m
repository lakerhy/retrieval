% ************************* Usage *****************************
% calculate the average snow temperature  
% ************************* Input *****************************
% intput_file of snow and ice with temp [K], density [kg/m3], thickness [cm], cor.length [mm],  salinity [per mil], ice/snow[1/0]
% snow_layernums: numbers of snow layers
% *************************Output*****************************

function[T_snow_FY,T_snow_MY,T_ice_FY,T_ice_MY]  = phytemp(p,snow_layernums)

  Tsea=271.35;      %[K] 
  ks=0.31;
  k0=2.034;         %[W/m/K]
  beta=0.117;        %[W/m/per mil] 

  Ts_snow = p(1,1); 
  roi_snow = p(4,1);
  di_ice_MY = 3; %m
  di_ice_FY = 1; %m
  di_ice = [di_ice_FY,di_ice_MY] % [m]
  di_snow = p(2,1)/100 %[m]
  pci_snow = p(3,1);
  sal_MY = 2.5;
  sal_FY = 7;
  sal = [sal_FY,sal_MY];
  
  % Range of validity of surface tempurature of snow surface
  if Ts_snow > 271.35
    Ts_snow = 271.35;
  elseif Ts_snow < 229.95
    Ts_snow = 229.95;
  end

  T(1) = Tsea;
  for j = 1:2 % FY 1 and MY 2
    
    for i = 1:80
      %    ki = k0 + (beta*sal_ice)/(T(1)+273.15); %conductivity of ice 
      ki = k0 + (beta*sal(j))/T(1); %conductivity of ice 
      slope = ki*di_ice(j)/ks/di_snow; %Ratio between the temperature gradients in ice and snow
      T0(1) = (Ts_snow + slope*Tsea)/(1 + slope); %Temperature of the interface of
                                                  %snow and ice
      slope2 = (T0(1)-Tsea)/di_ice(j); % temperature gradient of ice layer
      if abs(slope2) > 0.1
        slope2 = slope2/abs(slope2)*0.1;
        T0(1) = Tsea+slope2*di_ice(j);
      end
      slope3 = (Ts_snow-T0(1))/di_snow;  % temperture slope of snow layer
      T(1) = Tsea+slope2*di_ice(j)/2;
      T(2) = T0(1)+slope3*di_snow/2;
      
      T(1) = (Tsea+T(1))/2;  % mean value of ice temperature
                             %T0(2)=T0(1)
      if (i>1 && abs(T0(1)-T0(2))<0.01) 
        break
      end
      T0(2)=T0(1);
    end 
    
    T=T0(1);
    T_ice=T(1);
    
    T0 = abs(T_ice - 273.15);
    if T0<1
      W_ice = 1e-3*sal(j)*49.717;
    else 
      W_ice = 1e-3*sal(j)*(-49.185/T_ice+0.532);
    end
    delta_T1 = (T - Ts_snow)/snow_layernums; 
    delta_T2 = Tsea - T;
    Ti(1) = Tsea - delta_T2/2;
    for i = 1:(snow_layernums+1)
      Ti(i+1) = T - (i-1)*delta_T1;
    end
    for i = 1:(snow_layernums)
      Ti(i+1) = (Ti(i+1) + Ti(i+2))/2;
    end
    T_av(j,:) = Ti;
  end
  T_snow_FY =  T_av(1,2);
  T_ice_FY  =  T_av(1,1);
  T_snow_MY =  T_av(2,2);
  T_ice_MY  =  T_av(2,1);

  
  

  

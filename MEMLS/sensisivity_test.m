function[deltaTb]= sensisivity_test(ice_type,x)

  if ice_type == 1 % MY
    input = load('input_MY');
  else  % type=2 FY
    input = load('input_FY');
  end
  V = 2;
  deltaV = 0.5;
  Ts_MY = input(2,1);
  deltaTs = 1;
  di_snow = input(2,3);
  deltaDi_snow = 5;
  roi_snow = input(2,2);
  deltaRoi = 1;
  pci_snow = input(2,4);
  deltaPci = 0.05;
  di_ice = input(1,3);
  deltaDi_ice =1;
  sal = input(1,5);
  deltaSal  = 1;
  
  p = [Ts_MY; V;di_snow;roi_snow;pci_snow;di_ice;sal ];        

  deltaTb = zeros(10,1);
  deltaE = zeros(10,1);
  freq = [6.9,10.7,18.7,23.8,36.5];
  % For Ts
  if x==1
    p2 = [Ts_MY+deltaTs; V;di_snow;roi_snow;pci_snow;di_ice;sal];
    %     deltaTb =
    %     abs(fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type)-fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type))/deltaTs;    
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type)
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type)
    deltaTb = abs(A1-A2)/deltaTs  
    deltaE = abs(B1-B2)/deltaTs 
    
    % Plot
    figure
    subplot(2,1,1)
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
    xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by surface temperature vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
    xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by surface temperature vs Freq');
    legend('e_v','e_h',2);
    hold off
    
    
  end
  % For V
  if x==2
    p2 = [Ts_MY; V+deltaV;di_snow;roi_snow;pci_snow;di_ice;sal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7, 1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2) /deltaV;
    deltaE = abs(B1-B2)/deltaV;
    % Plot without the MEMLS model
    if ice_type==1
    p1a = [0, V, 0, 0, Ts_MY, 1,1];
    p2a = [0, V+deltaV, 0, 0, Ts_MY, 1,1];
    else 
      
    p1a = [0, V, 0, 0, Ts_MY, 1,0];
    p2a = [0, V+deltaV, 0, 0, Ts_MY, 1,0];
    
    end    
    
    deltaTba = abs(FW_funktion2_is(p1a(1,1),p1a(1,2),p1a(1,3),p1a(1,4),p1a(1,5),p1a(1,6),p1a(1, ...
                                                      7))- ...
                   FW_funktion2_is(p2a(1,1),p2a(1,2),p2a(1,3),p2a(1,4),p2a(1, ...
                                                      5),p2a(1,6),p2a(1,7)) ) /deltaV;
    figure
%    subplot(2,1,1)
    grid on
    hold on
    %V channel
    
    plot(freq,deltaTba(1:2:10),'--mo','LineWidth',1.5);
    %H channel
    plot(freq,deltaTba(2:2:10),'--go','LineWidth',1.5);
    xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by water vapor content vs Freq with/out MEMLS');
    % legend('Tbv_o','Tbh_o',2);
    hold off
    % Plot
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);

    legend('Tbv_o','Tbh_o','Tbv','Tbh',4);
    hold off
     
%     subplot(2,1,2)
%     grid on
%     hold on
%     %V channel
%     plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
%     %H channel
%     plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
%     xlabel('Frequency [GHz]');
%     title('FY: Sensitivity of Emissivity by vapor content  vs Freq');
%     legend('e_v','e_h',2);
%     hold off

  end
  % For di_snow
  if x==3
    p2 = [Ts_MY; V;di_snow+deltaDi_snow;roi_snow;pci_snow;di_ice;sal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2)/deltaDi_snow;    
    deltaE = abs(B1-B2)/deltaDi_snow;    

    % Plot
    figure
    subplot(2,1,1)
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
    xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by snow depth vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
    xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by snow depth  vs Freq');
    legend('e_v','e_h',2);
    hold off
  end
  % For roi_snow
  
  if x==4
    p2 = [Ts_MY; V;di_snow;roi_snow+deltaRoi;pci_snow;di_ice;sal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2)/deltaRoi;
    deltaE = abs(B1-B2)/deltaRoi;
    % Plot
    figure
    subplot(2,1,1)
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by snow density vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    
    
       subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by  snow density vs Freq');
    legend('E_v','E_h',2);
    hold off

  end
  % For pci_snow
  if x==5
    p2 = [Ts_MY; V;di_snow;roi_snow;pci_snow+deltaPci;di_ice;sal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2)/deltaPci;
    deltaE = abs(B1-B2)/deltaPci;
    % Plot
    figure
        subplot(2,1,1)
    grid on
    hold on
    %V channel
    
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by correlation length of snow vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    
    
        subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by correlation length of snow vs Freq');
    legend('e_v','e_h',2);
    hold off

  end
  % For di_ice
  if x==6
    p2 = [Ts_MY; V;di_snow;roi_snow;pci_snow;di_ice+deltaDi_ice;sal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2)/deltaDi_ice;
    deltaE = abs(B1-B2)/deltaDi_ice;
    % Plot
    figure
        subplot(2,1,1)
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');

    title('FY: Sensitivity of Tb by thickness of ice vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    
    
        subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by thickness of ice  vs Freq');
    legend('e_v','e_h',2);
    hold off

  end
  % For sal
  if x==7
    p2 = [Ts_MY; V;di_snow;roi_snow;pci_snow;di_ice;sal+deltaSal ];
    [A1,B1] = fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),ice_type);
    [A2,B2] = fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),ice_type);
    deltaTb = abs(A1-A2)/deltaSal;
    deltaE = abs(B1-B2)/deltaSal;
    
    
    
    % Plot
    figure
        subplot(2,1,1)
    grid on
    hold on
    %V channel
    plot(freq,deltaTb(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaTb(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Tb by salinity vs Freq');
    legend('Tbv','Tbh',2);
    hold off
    
        subplot(2,1,2)
    grid on
    hold on
    %V channel
    plot(freq,deltaE(1:2:10),'-rs','LineWidth',1.5);
    %H channel
    plot(freq,deltaE(2:2:10),'-bs','LineWidth',1.5);
xlabel('Frequency [GHz]');
    title('FY: Sensitivity of Emissivity by salinity vs Freq');
    legend('e_v','e_h',2);
    hold off

  end

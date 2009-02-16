% type 1/2 MY/FY
function Tb_plot_fun(a,type,p,output)
 % V=0.2;

   freq = [6.9,10.7,18.7,23.8,36.5];

  if type==1 % MY
             %     discretize('input_MY',snow_layernums,'syntes.disc.MY',1);
             %     Tb_MY = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
             %                     'syntes.disc.MY',4);
             %     emissivity_MY = emi('syntes.disc.MY',4)
     Ts_MY = p(1,1);
        V = p(2,1)
    di_snow = p(3,1);
    roi_snow = p(4,1);
    pci_snow = p(5,1);
    di_ice = p(6,1);
    sal = p(7,1);
    p_input = [Ts_MY; V;di_snow;roi_snow;pci_snow;di_ice;sal ];    
    TB_MY=fw_fun2(p_input(1,1),p_input(2,1),p_input(3,1),p_input(4,1),p_input(5,1),p_input(6,1),p_input(7,1),type);
    %% ################### PLOT ############################### 
 subplot(2,2,a)
    grid on
    hold on
  
    plot(freq,TB_MY(1:2:10),'-rs');
    plot(freq,TB_MY(2:2:10),'-bs');
    plot(freq,output(1:2:10),'--ms');%V channel
    plot(freq,output(2:2:10),'--gs');%H channe
  %  xlabel('Frequency[GHz]');
   % ylabel('Brightness Temperature[K]');
%    title('MY: Brightness temperature vs. Frequency in Jan of area 33 ');
    legend('Tbv','Tbh','Tbv_a','Tbh_a',2);
    
  end
  
  if type==2 % FY 
         Ts_FY = p(1,1);
     V = p(2,1)
    di_snow = p(3,1);
    roi_snow = p(4,1);
    pci_snow = p(5,1);
    di_ice = p(6,1);
    sal = p(7,1);
    p_input = [Ts_FY; V;di_snow;roi_snow;pci_snow;di_ice;sal ];    
    TB_FY=fw_fun2(p_input(1,1),p_input(2,1),p_input(3,1),p_input(4,1),p_input(5,1),p_input(6,1),p_input(7,1),type);
    %% ################### PLOT ############################### 
 subplot(2,2,a)
    grid on
    hold on
    
    plot(freq,TB_FY(1:2:10),'-rs');
    plot(freq,TB_FY(2:2:10),'-bs');
    plot(freq,output(1:2:10),'--ms');%V channel
    plot(freq,output(2:2:10),'--gs');%H channe
    xlabel('Frequency [GHz]');
    ylabel('Brightness Temperature [K]');
    %title('FY: Brightness temperature vs. Frequency in Jan of area 58 ');
    legend('Tbv','Tbh','Tbv_a','Tbh_a',2);

  end


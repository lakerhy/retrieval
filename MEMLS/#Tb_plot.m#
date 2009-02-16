% type 1/2 MY/FY
function Tb_plot(V,type,par)
  
  output = zeros(1,10);
  TB_MY= zeros(1,10);
  TB_FY=zeros(1,10);
  

  
  % V=0.2;
  
  freq = [6.9,10.7,18.7,23.8,36.5];
  if type==1 % MY
             %     discretize('input_MY',snow_layernums,'syntes.disc.MY',1);
             %     Tb_MY = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
             %                     'syntes.disc.MY',4);
             %     emissivity_MY = emi('syntes.disc.MY',4)
    cd Data;
    output = average_T('area.33.sort',1);
    cd ..;
    TB_MY=fw_fun2(par,type,V);
 
    %% ################### PLOT ############################### 
    figure
    grid on
    hold on
    plot(freq,TB_MY(1:2:10),'-ro', 'LineWidth',1.5);
    plot(freq,TB_MY(2:2:10),'-bo', 'LineWidth',1.5);
    plot(freq,output(1:2:10),'--ms', 'LineWidth',1.5);%V channel
    plot(freq,output(2:2:10),'--gs', 'LineWidth',1.5);%H channe
    xlabel('Frequency [GHz]');
    ylabel('Brightness Temperature [K]');
    title('MY: Brightness temperature vs. Frequency in Jan of area 33 ');
    legend('Tbv','Tbh','Tbv_a','Tbh_a',2);
    hold off
  end
  
  if type==2 % FY 
    cd Data;
    output = average_T('area.58s.sort',2);
    cd ..;
    
    TB_FY=fw_fun2(par,type,V);
    %% ################### PLOT ############################### 
    figure
    grid on
    hold on
    
    plot(freq,TB_FY(1:2:10),'-ro', 'LineWidth',1.5);
    plot(freq,TB_FY(2:2:10),'-bo', 'LineWidth',1.5);
    plot(freq,output(1:2:10),'--ms', 'LineWidth',1.5);%V channel
    plot(freq,output(2:2:10),'--gs', 'LineWidth',1.5);%H channe
    xlabel('Frequency [GHz]');
    ylabel('Brightness Temperature [K]');
    title('FY: Brightness temperature vs. Frequency in Jan of area 58 ');
    legend('Tbv','Tbh','Tbv_a','Tbh_a',2);
    hold off
  end

  

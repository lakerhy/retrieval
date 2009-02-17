function[] = memlsamsr (type)

  frequency = [6.9,10.7,18.7,23.8,36.5]; 
  
  if type ==1 %FY
%     cd ../MEMLS
%     %FY=load('FY.profile.1');
%     FY=load('FY.profile');
%     Tb_FY=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', FY);    
%     Tb_FY_V=Tb_FY(:,1);
%     Tb_FY_H=Tb_FY(:,2);
    
           %FY=load('FY.profile.1');
           cd ../MEMLS
    FY=load('FY.profile.1');
    cd ../tune
    Tb_FY=fwtune(FY)
    Tb_FY_V=Tb_FY(:,1)
    Tb_FY_H=Tb_FY(:,2)

    cd ../inversion
    Tb_amsr = amsr('area.58.2005.sort');
    Tb_amsr_V = Tb_amsr(1:2:10);
    Tb_amsr_H = Tb_amsr(2:2:10);
    Tb_amsr2 = amsr('area.46.2005.sort');
    Tb_amsr_V2 = Tb_amsr2(1:2:10);
    Tb_amsr_H2 = Tb_amsr2(2:2:10);
    cd ../tune
    
    figure
    
    subplot(2,1,1)
    hold on
    grid on
    plot(frequency,Tb_FY_V,'r+');
    plot(frequency,Tb_amsr_V,'--g');
    plot(frequency,Tb_amsr_V2,'--b');
    legend('sim','amsr58','smsr46');
    ylabel('Tb_V');
    xlabel('frequency');
    title('AMSR and simulated.');

    subplot(2,1,2)
    hold on
    grid on
    plot(frequency,Tb_FY_H,'r+');
    plot(frequency,Tb_amsr_H,'g--');
    plot(frequency,Tb_amsr_H2,'--b');
    legend('sim','amsr58','smsr46');
    ylabel('Tb_H');
    xlabel('frequency');
    title('AMSR and simulated');
  end

  if type ==2 %MY
    cd ../MEMLS
    %    MY=load('MY.profile.2');
    MY=load('MY.profile.1');
 %    [Tb_MY]=icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', MY)
%     T_MY_V=Tb_MY(:,1);
%     T_MY_H=Tb_MY(:,2);
               cd ../tune
     Tb_MY=fwtune(MY)
    T_MY_V=Tb_MY(:,1)
    T_MY_H=Tb_MY(:,2)

    
    
    cd ../inversion
    Tb_amsr = amsr('area.33.2005.sort');
    Tb_amsr_V = Tb_amsr(1:2:10);
    Tb_amsr_H = Tb_amsr(2:2:10);
    Tb_amsr2 = amsr('area.32.2005.sort');
    Tb_amsr_V2 = Tb_amsr2(1:2:10);
    Tb_amsr_H2 = Tb_amsr2(2:2:10);
    cd ../tune
    
    figure
    
    subplot(2,1,1)
    hold on
    grid on
    plot(frequency,T_MY_V,'r+');
    plot(frequency,Tb_amsr_V,'--g');
    plot(frequency,Tb_amsr_V2,'--b');
    legend('sim','amsr33','amsr32');
    ylabel('Tb_V');
    xlabel('frequency');
    title('AMSR and simulated');

    subplot(2,1,2)
    hold on
    grid on
    plot(frequency,T_MY_H,'r+');
    plot(frequency,Tb_amsr_H,'--g');
    plot(frequency,Tb_amsr_H2,'--b');
    legend('sim','amsr33','amsr32');
    ylabel('Tb_H');
    xlabel('frequency');
    title('AMSR and simulated');            
    
  end


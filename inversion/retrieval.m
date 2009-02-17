function[] = retrieval(type)

  if type == 1 
    % Inversion evaluation
    %p=[250,18,0.1,320,0]'; 
    p=[231,30,0.14,350,0]';
    Tb = fw(p);
    
  end

  if type == 2
    % FY data 58s    
 %   [Tb_AMSR]=amsr('area.46.2005.sort');
      [Tb_AMSR]=amsr('area.58.2005.sort');
    Tb = Tb_AMSR;
  end
  
  if type ==3
    % MY data 33s    
%    [Tb_AMSR]=amsr('area.33.2005.sort');
       [Tb_AMSR]=amsr('area.33.2005.sort');
    Tb = Tb_AMSR;
  end

  [p_est,S_std,Sp_std]=inversion(Tb);
   Tb_est = fw(p_est);
   %[Tb_AMSR]=amsr('area.33.2005.sort');
    Tb_amsr_V = Tb(1:2:10);
    Tb_amsr_H = Tb(2:2:10);
            
    figure
    frequency = [6.9,10.7,18.7,23.8,36.5]; 
    subplot(2,1,1)
    hold on
    grid on
    plot(frequency,Tb_est(1:2:10),'r+');
    plot(frequency,Tb_amsr_V,'--g');
   %plot(frequency,Tb_AMSR(1:2:10),'--b');
    %%%%%legend('est','amsr33');
    ylabel('Tb_V');
    xlabel('frequency');
    title('AMSR and estimation.');

    subplot(2,1,2)
    hold on
    grid on
    plot(frequency,Tb_est(2:2:10),'r+');
    plot(frequency,Tb_amsr_H,'g--');
    %   plot(frequency,Tb_AMSR(2:2:10),'--b')
    %legend('est','amsr33');
    ylabel('Tb_H');
    xlabel('frequency');
    title('AMSR and estimation');
p_est;
S_std;
Sp_std;

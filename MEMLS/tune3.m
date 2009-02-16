
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             K               kg/m^-3 cm 
%%      num,   T          W   ,  roi   ,    di ,     pci        sal  type
%%  par = [1 , Tsea ,   0.0026, roi_ice,   di_ice,   0.35,      7.5,  1;
%%         2 , Ts_snow ,0.0000, roi_snow2, di_snow2, pci_snow2, 0,    0;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = tune3(n)
% Introducing the slushy layer or not
% Tb under different profiles on same plot  
% n: 1/2/3 snow density/correlation length/depth

  par_in = load('input_FY') ;
  par_in_S = load('input_FY_S');
  density_step = 40;
  delta_density = linspace(100,500,density_step);
  pci_step = 40;
  delta_pci = linspace(0.01,0.2,pci_step);
  di_step = 50;
  delta_di = linspace(5,55,di_step);

  output = zeros(1,10);
  TB_MY= zeros(1,10);
  TB_FY=zeros(1,10);
  freq = [6.9,10.7,18.7,23.8,36.5];
  cd Data;
  output = average_T('area.58s.sort',2);
  cd ..;

  if n==1
    TB_FY_6G = zeros(2,density_step);
    TB_FY_11G = zeros(2,density_step);
    TB_FY_19G = zeros(2,density_step);
    TB_FY_24G = zeros(2,density_step);
    TB_FY_37G = zeros(2,density_step);
    TB_FY_6G_S = zeros(2,density_step);
    TB_FY_11G_S = zeros(2,density_step);
    TB_FY_19G_S = zeros(2,density_step);
    TB_FY_24G_S = zeros(2,density_step);
    TB_FY_37G_S = zeros(2,density_step);
    for i = 1:density_step
      par_in(3,4) = delta_density(i);
      par_in_S(3,4) = delta_density(i);
      par_out =  homo(par_in,10,2);
      par_out = par_out';
      par_out_S =  homo(par_in_S,10,2);
      par_out_S = par_out_S';
      TB_FY=fw_fun2(par_out,2,0.2);
      TB_FY_S=fw_fun2(par_out_S,2,0.2);
      TB_FY_6G(:,i) = [TB_FY(1,1);TB_FY(1,2)];% V;H
      TB_FY_11G(:,i) = [TB_FY(1,3);TB_FY(1,4)];% V;H
      TB_FY_19G(:,i) = [TB_FY(1,5);TB_FY(1,6)];% V;H
      TB_FY_24G(:,i) = [TB_FY(1,7);TB_FY(1,8)];% V;H
      TB_FY_37G(:,i) = [TB_FY(1,9);TB_FY(1,10)];% V;H
      TB_FY_6G_S(:,i) = [TB_FY_S(1,1);TB_FY_S(1,2)];% V;H
      TB_FY_11G_S(:,i) = [TB_FY_S(1,3);TB_FY_S(1,4)];% V;H
      TB_FY_19G_S(:,i) = [TB_FY_S(1,5);TB_FY_S(1,6)];% V;H
      TB_FY_24G_S(:,i) = [TB_FY_S(1,7);TB_FY_S(1,8)];% V;H
      TB_FY_37G_S(:,i) = [TB_FY_S(1,9);TB_FY_S(1,10)];% V;H
    end
    
    %%%%%%%%%V channel Snow Density%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on;
    plot(delta_density,abs(TB_FY_6G(1,:)- output(1,1)  ),'r');
    plot(delta_density,abs(TB_FY_11G(1,:)-output(1,3)  ),'g');
    plot(delta_density,abs(TB_FY_19G(1,:)-output(1,5)  ),'b');
    plot(delta_density,abs(TB_FY_24G(1,:)-output(1,7)  ),'c');
    plot(delta_density,abs(TB_FY_37G(1,:)-output(1,9)  ),'m');
    plot(delta_density,abs(TB_FY_6G_S(1,:)-output(1,1) ),'--r');
    plot(delta_density,abs(TB_FY_11G_S(1,:)-output(1,3)),'--g');
    plot(delta_density,abs(TB_FY_19G_S(1,:)-output(1,5)),'--b');
    plot(delta_density,abs(TB_FY_24G_S(1,:)-output(1,7)),'--c');
    plot(delta_density,abs(TB_FY_37G_S(1,:)-output(1,9)),'--m');
    legend('6V','11V','19V','24V','37V');
    title('TB_V difference vs snow density with/out slushy layer');
    xlabel('Snow density [kg/m^{3}]');
    ylabel('|TB_{memls}_V-TB_V| [K]');
    hold off
    
    %%%%%%%%%H channel Snow Density%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_density,abs(TB_FY_6G(2,:)-output(1,2)   ),'r');
    plot(delta_density,abs(TB_FY_11G(2,:)-output(1,4)  ) ,'g');    
    plot(delta_density,abs(TB_FY_19G(2,:)-output(1,6)  ) ,'b');
    plot(delta_density,abs(TB_FY_24G(2,:)-output(1,8) ) ,'c');
    plot(delta_density,abs(TB_FY_37G(2,:)-output(1,10) ) ,'m');
    plot(delta_density,abs(TB_FY_6G_S(2,:)-output(1,2) ) ,'--r');
    plot(delta_density,abs(TB_FY_11G_S(2,:)-output(1,4)),'--g');
    plot(delta_density,abs(TB_FY_19G_S(2,:)-output(1,6)),'--b');
    plot(delta_density,abs(TB_FY_24G_S(2,:)-output(1,8)),'--c');
    plot(delta_density,abs(TB_FY_37G_S(2,:)-output(1,10)),'--m');
    %legend('6GH','11GH','19GH','24GH','37GH','6GH_S','11GH_S','19GH_S', ...
    %       '24GH_S','37GH_S');
legend('6H','11H','19H','24H','37H');
    title('TB_H difference vs snow density with/out slushy layer');
    xlabel('Snow density [kg/m^{3}]');
    ylabel('|TB_{memls}_H-TB_H| [K]');
    hold off
    
    %%%%%%%%%V-H polarization Snow Density%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_density,abs((TB_FY_6G(1,:)-TB_FY_6G(2,:))-(output(1,1)-output(1,2))),'r');
    plot(delta_density,abs((TB_FY_11G(1,:)-TB_FY_11G(2,:))-(output(1,3)-output(1,4))),'g');
    plot(delta_density,abs((TB_FY_19G(1,:)-TB_FY_19G(2,:))-(output(1,5)-output(1,6))),'b');
    plot(delta_density,abs((TB_FY_24G(1,:)-TB_FY_24G(2,:))-(output(1,7)-output(1,8))),'c');
    plot(delta_density,abs((TB_FY_37G(1,:)-TB_FY_37G(2,:))-(output(1,9)-output(1,10))),'m');
    plot(delta_density,abs(( TB_FY_6G_S(1,:) -TB_FY_6G_S(2,:))-(output(1,1)-output(1,2))),'--r');
    plot(delta_density,abs((TB_FY_11G_S(1,:)-TB_FY_11G_S(2,:))-(output(1,3)-output(1,4))),'--g');
    plot(delta_density,abs((TB_FY_19G_S(1,:)-TB_FY_19G_S(2,:))-(output(1,5)-output(1,6))),'--b');
    plot(delta_density,abs((TB_FY_24G_S(1,:)-TB_FY_24G_S(2,:))-(output(1,7)-output(1,8))),'--c');
    plot(delta_density,abs((TB_FY_37G_S(1,:)-TB_FY_37G_S(2,:))-(output(1,9)-output(1,10))),'--m');
    %legend('6GV-H','11GV-H','19GV-H','24GV-H','37GV-H','6GV-H_S',['11GV-' ...
     %                   'H_S'],'19GV-H_S','24GV-H_S','37GV-H_S');

    legend('6V-H','11V-H','19V-H','24V-H','37V-H')
    title('Polarization difference vs snow density with/out slushy layer');
    xlabel('Snow density [kg/m^{3}]');
    ylabel('|TB_{memls}_{V-H}-TB_{V-H}| [K]');
    hold off
  end
  
  if n==2
    TB_FY_6G = zeros(2, pci_step);
    TB_FY_11G = zeros(2,pci_step);
    TB_FY_19G = zeros(2,pci_step);
    TB_FY_24G = zeros(2,pci_step);
    TB_FY_37G = zeros(2,pci_step);
    TB_FY_6G_S = zeros(2, pci_step);
    TB_FY_11G_S = zeros(2,pci_step);
    TB_FY_19G_S = zeros(2,pci_step);
    TB_FY_24G_S = zeros(2,pci_step);
    TB_FY_37G_S = zeros(2,pci_step);
    
    for i = 1:pci_step
      par_in(3,6) = delta_pci(i);
      par_in_S(3,6) = delta_pci(i);
      par_out =  homo(par_in,10,2);
      par_out = par_out';
      par_out_S =  homo(par_in_S,10,2);
      par_out_S = par_out_S';
      TB_FY=fw_fun2(par_out,2,0.2);
      TB_FY_S=fw_fun2(par_out_S,2,0.2);
      TB_FY_6G(:,i) = [TB_FY(1,1);TB_FY(1,2)];% V;H
      TB_FY_11G(:,i) = [TB_FY(1,3);TB_FY(1,4)];% V;H
      TB_FY_19G(:,i) = [TB_FY(1,5);TB_FY(1,6)];% V;H
      TB_FY_24G(:,i) = [TB_FY(1,7);TB_FY(1,8)];% V;H
      TB_FY_37G(:,i) = [TB_FY(1,9);TB_FY(1,10)];% V;H
      TB_FY_6G_S(:,i) = [TB_FY_S(1,1);TB_FY_S(1,2)];% V;H
      TB_FY_11G_S(:,i) = [TB_FY_S(1,3);TB_FY_S(1,4)];% V;H
      TB_FY_19G_S(:,i) = [TB_FY_S(1,5);TB_FY_S(1,6)];% V;H
      TB_FY_24G_S(:,i) = [TB_FY_S(1,7);TB_FY_S(1,8)];% V;H
      TB_FY_37G_S(:,i) = [TB_FY_S(1,9);TB_FY_S(1,10)];% V;H
    end
    
    %%%%%%%%%%%%%%%%%%%%%V channel Snow Correlation length%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on;
    plot(delta_pci,abs(TB_FY_6G(1,:)-output(1,1)  ), 'r');
    plot(delta_pci,abs(TB_FY_11G(1,:)-output(1,3) ),'g');
    plot(delta_pci,abs(TB_FY_19G(1,:)-output(1,5) ),'b');
    plot(delta_pci,abs(TB_FY_24G(1,:)-output(1,7) ),'c');
    plot(delta_pci,abs(TB_FY_37G(1,:)-output(1,9) ),'m');
    plot(delta_pci,abs(TB_FY_6G_S(1,:)-output(1,1 )),'--r');
    plot(delta_pci,abs(TB_FY_11G_S(1,:)-output(1,3)),'--g');
    plot(delta_pci,abs(TB_FY_19G_S(1,:)-output(1,5)),'--b');
    plot(delta_pci,abs(TB_FY_24G_S(1,:)-output(1,7)),'--c');
    plot(delta_pci,abs(TB_FY_37G_S(1,:)-output(1,9)),'--m');
legend('6V','11V','19V','24V','37V');
    title('TB_V difference vs correlation length with/out slushy layer');
    xlabel('snow correlation length [mm]');
    ylabel('|TB_{memls}_V-TB_V| [K]');
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%H channel Snow Correlation length%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_pci,abs(TB_FY_6G(2,:)-output(1,2)  ),'r');
    plot(delta_pci,abs(TB_FY_11G(2,:)-output(1,4) ),'g');    
    plot(delta_pci,abs(TB_FY_19G(2,:)-output(1,6) ),'b');
    plot(delta_pci,abs(TB_FY_24G(2,:)-output(1,8)),'c');
    plot(delta_pci,abs(TB_FY_37G(2,:)-output(1,10)),'m');
    plot(delta_pci,abs(TB_FY_6G_S(2,:)-output(1,2)),'--r');
    plot(delta_pci,abs(TB_FY_11G_S(2,:)-output(1,4)),'--g');
    plot(delta_pci,abs(TB_FY_19G_S(2,:)-output(1,6)),'--b');
    plot(delta_pci,abs(TB_FY_24G_S(2,:)-output(1,8)),'--c');
    plot(delta_pci,abs(TB_FY_37G_S(2,:)-output(1,10)),'--m');
   legend('6H','11H','19H','24H','37V');
   title('TB_H vs correlation length with/out slushy layer');
    xlabel('snow correlation length [mm]');
    
    ylabel('|TB_{memls}_H-TB_H| [K]');
    hold off

    %%%%%%%%%%%%%%%%%%%%%V-H polariztion Snow Correlation length%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_pci,abs((TB_FY_6G(1,:)-TB_FY_6G(2,:))-(output(1,1)-output(1,2))),'r');
    plot(delta_pci,abs((TB_FY_11G(1,:)-TB_FY_11G(2,:))-(output(1,3)-output(1,4))),'g');
    plot(delta_pci,abs((TB_FY_19G(1,:)-TB_FY_19G(2,:))-(output(1,5)-output(1,6))),'b');
    plot(delta_pci,abs((TB_FY_24G(1,:)-TB_FY_24G(2,:))-(output(1,7)-output(1,8))),'c');
    plot(delta_pci,abs((TB_FY_37G(1,:)-TB_FY_37G(2,:))-(output(1,9)-output(1,10))),'m');
    plot(delta_pci,abs((TB_FY_6G_S(1,:) -TB_FY_6G_S(2,:))-(output(1,1)-output(1,2))),'--r');
    plot(delta_pci,abs((TB_FY_11G_S(1,:)-TB_FY_11G_S(2,:))-(output(1,3)-output(1,4))),'--g');
    plot(delta_pci,abs((TB_FY_19G_S(1,:)-TB_FY_19G_S(2,:))-(output(1,5)-output(1,6))),'--b');
    plot(delta_pci,abs((TB_FY_24G_S(1,:)-TB_FY_24G_S(2,:))-(output(1,7)-output(1,8))),'--c');
    plot(delta_pci,abs((TB_FY_37G_S(1,:)-TB_FY_37G_S(2,:))-(output(1,9)-output(1,10))),'--m');
   legend('6V-H','11V-H','19V-H','24V-H','37V-H');
       title('Polarization difference vs snow correlation length with/out slushy layer');
    xlabel('Snow correlation length [mm]');
    ylabel('|TB_{memls}_{V-H}-TB_{V-H}| [K]');
    hold off
  end

  if n==3
    TB_FY_6G = zeros(2, di_step);
    TB_FY_11G = zeros(2,di_step);
    TB_FY_19G = zeros(2,di_step);
    TB_FY_24G = zeros(2,di_step);
    TB_FY_37G = zeros(2,di_step);
    TB_FY_6G_S = zeros(2,di_step);
    TB_FY_11G_S = zeros(2,di_step);
    TB_FY_19G_S = zeros(2,di_step);
    TB_FY_24G_S = zeros(2,di_step);
    TB_FY_37G_S = zeros(2,di_step);
    for i = 1:di_step
      par_in(3,5) = delta_di(i);
      par_in_S(3,5) = delta_di(i);
      par_out =  homo(par_in,10,2);
      par_out = par_out';
      par_out_S =  homo(par_in_S,10,2);
      par_out_S = par_out_S';
      TB_FY=fw_fun2(par_out,2,0.2);
      TB_FY_S=fw_fun2(par_out_S,2,0.2);
      TB_FY_6G_S(:,i) = [TB_FY_S(1,1);TB_FY_S(1,2)];% V;H
      TB_FY_11G_S(:,i) = [TB_FY_S(1,3);TB_FY_S(1,4)];% V;H
      TB_FY_19G_S(:,i) = [TB_FY_S(1,5);TB_FY_S(1,6)];% V;H
      TB_FY_24G_S(:,i) = [TB_FY_S(1,7);TB_FY_S(1,8)];% V;H
      TB_FY_37G_S(:,i) = [TB_FY_S(1,9);TB_FY_S(1,10)];% V;H
      TB_FY_6G(:,i) = [TB_FY(1,1);TB_FY(1,2)];% V;H
      TB_FY_11G(:,i) = [TB_FY(1,3);TB_FY(1,4)];% V;H
      TB_FY_19G(:,i) = [TB_FY(1,5);TB_FY(1,6)];% V;H
      TB_FY_24G(:,i) = [TB_FY(1,7);TB_FY(1,8)];% V;H
      TB_FY_37G(:,i) = [TB_FY(1,9);TB_FY(1,10)];% V;H
    end
    
    %%%%%%%%%%%%%%%%V-channel di%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on;
    
    plot(delta_di,abs(TB_FY_6G(1,:)-output(1,1)  ),'r');
    plot(delta_di,abs(TB_FY_11G(1,:)-output(1,3) ),'g');
    plot(delta_di,abs(TB_FY_19G(1,:)-output(1,5) ),'b');
    plot(delta_di,abs(TB_FY_24G(1,:)-output(1,7) ),'c');
    plot(delta_di,abs(TB_FY_37G(1,:)-output(1,9)),'m');
    plot(delta_di,abs(TB_FY_6G_S(1,:)-output(1,1)),'--r');
    plot(delta_di,abs(TB_FY_11G_S(1,:)-output(1,3)),'--g');
    plot(delta_di,abs(TB_FY_19G_S(1,:)-output(1,5)),'--b');
    plot(delta_di,abs(TB_FY_24G_S(1,:)-output(1,7)),'--c');
    plot(delta_di,abs(TB_FY_37G_S(1,:)-output(1,9)),'--m');
    legend('6V','11V','19V','24V','37V');
    title('TB_V difference vs snow depth with/out slushy layer');
    xlabel('Snow depth [cm]');
    ylabel('|TB_{memls}_V-TB_V| [K]');
    hold off
    
    %%%%%%%%%%%%%%%%%%%H-channel di%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_di,abs(TB_FY_6G(2,:)-output(1,2)   ),'r');
    plot(delta_di,abs(TB_FY_11G(2,:)-output(1,4)  ),'g');    
    plot(delta_di,abs(TB_FY_19G(2,:)-output(1,6)  ), 'b');
    plot(delta_di,abs(TB_FY_24G(2,:)-output(1,8) ),'c');
    plot(delta_di,abs(TB_FY_37G(2,:)-output(1,10) ),'m');
    plot(delta_di,abs(TB_FY_6G_S(2,:)-output(1,2) ),'--r');
    plot(delta_di,abs(TB_FY_11G_S(2,:)-output(1,4 )),'--g');
    plot(delta_di,abs(TB_FY_19G_S(2,:)-output(1,6 )),'--b');
    plot(delta_di,abs(TB_FY_24G_S(2,:)-output(1,8 )),'--c');
    plot(delta_di,abs(TB_FY_37G_S(2,:)-output(1,10)),'--m');
    title('TB_H difference vs snow density with/out slushy layer');
    legend('6H','11H','19H','24H','37H');
    
          xlabel('Snow depth [cm]');
    ylabel('|TB_{memls}_H-TB_H| [K]');
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%V-H polarization di%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    grid on
    hold on
    plot(delta_di,abs((TB_FY_6G(1,:)-TB_FY_6G(2,:))-(output(1,1)-output(1,2))),'r');
    plot(delta_di,abs((TB_FY_11G(1,:)-TB_FY_11G(2,:))-(output(1,3)-output(1,4))),'g');
    plot(delta_di,abs((TB_FY_19G(1,:)-TB_FY_19G(2,:))-(output(1,5)-output(1,6))),'b');
    plot(delta_di,abs((TB_FY_24G(1,:)-TB_FY_24G(2,:))-(output(1,7)-output(1,8))),'c');
    plot(delta_di,abs((TB_FY_37G(1,:)-TB_FY_37G(2,:))-(output(1,9)-output(1,10))),'m');
    plot(delta_di,abs(( TB_FY_6G_S(1,:) -TB_FY_6G_S(2,:))-(output(1,1)-output(1,2))),'--r');
    plot(delta_di,abs((TB_FY_11G_S(1,:)-TB_FY_11G_S(2,:))-(output(1,3)-output(1,4))),'--g');
    plot(delta_di,abs((TB_FY_19G_S(1,:)-TB_FY_19G_S(2,:))-(output(1,5)-output(1,6))),'--b');
    plot(delta_di,abs((TB_FY_24G_S(1,:)-TB_FY_24G_S(2,:))-(output(1,7)-output(1,8))),'--c');
    plot(delta_di,abs((TB_FY_37G_S(1,:)-TB_FY_37G_S(2,:))-(output(1,9)-output(1,10))),'--m');
    legend('6V-H','11V-H','19V-H','24V-H','37V-H');
    title('Polarization difference vs snow depth  length with/out slushy layer');
    xlabel('Snow depth [cm]');
    ylabel('|TB_{memls}_{V-H}-TB_{V-H}| [K]');
    hold off

  end

  

  
  
  
  
  
  
  
  
  
  
  
  

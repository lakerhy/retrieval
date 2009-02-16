%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             K               kg/m^-3 cm 
%%      num,   T          W   ,  roi   ,    di ,     pci        sal  type
%%  par = [1 , Tsea ,   0.0026, roi_ice,   di_ice,   0.35,      7.5,  1;
%%         2 , Ts_snow ,0.0000, roi_snow2, di_snow2, pci_snow2, 0,    0;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[hs,hs_AMSR] = tune4(n)
% Retrieval method using hs=2.9-782XGR_ice Markus eq(4)

  par_in = load('input_FY') ;
    par_in_S = load('input_FY_S') ;


  density_step = 35;
  delta_density = linspace(150,500,density_step);
  pci_step = 30;
  delta_pci = linspace(0.01,0.15,pci_step);
  di_step = 30;
  delta_di = linspace(10,70,di_step);

  output = zeros(1,10);
  TB_MY= zeros(1,10);
  TB_FY=zeros(1,10);


  freq = [6.9,10.7,18.7,23.8,36.5];

  cd Data;
  output = average_T('area.58s.sort',2);
  cd ..;

    
    if n==1
      
          
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
    
    GR = (TB_FY_37G(1,:)-TB_FY_19G(1,:))./(TB_FY_37G(1,:)+TB_FY_19G(1,:));    
    GR_S = (TB_FY_37G_S(1,:)-TB_FY_19G_S(1,:))./(TB_FY_37G_S(1,:)+ ...
                                                TB_FY_19G_S(1,:));
GR_AMSR = (output(1,9)-output(1,5))/(output(1,9)+output(1,5));


    hs= 2.9-782*GR;
    hs_S=2.9-782*GR_S;
    hs_AMSR= 2.9-782*GR_AMSR;
    
    figure
    grid on
    hold on;
    
    plot(delta_di,hs,'r');
    plot(delta_di,hs_S,'g');
    plot(delta_di,delta_di,'k');
    axis([10,70,10,35]);
%    plot(hs_AMSR*ones(1,30),delta_di,'k');
    title('Retrieval vs real snow depth with/out slushy layer');
    xlabel('Real snow depth [cm]');
    ylabel('Retrieval snow depth [cm]');
    legend('without','with');
    hold off 
    end
    
    
    if n==2
      
          
    TB_FY_6G = zeros(2, density_step);
    TB_FY_11G = zeros(2,density_step);
    TB_FY_19G = zeros(2,density_step);
    TB_FY_24G = zeros(2,density_step);
    TB_FY_37G = zeros(2,density_step);
    TB_FY_6G_S = zeros(2, density_step);
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
    
    GR = (TB_FY_37G(1,:)-TB_FY_19G(1,:))./(TB_FY_37G(1,:)+TB_FY_19G(1,:));    
    GR_S = (TB_FY_37G_S(1,:)-TB_FY_19G_S(1,:))./(TB_FY_37G_S(1,:)+ ...
                                                TB_FY_19G_S(1,:));
GR_AMSR = (output(1,9)-output(1,5))/(output(1,9)+output(1,5));


    hs= 2.9-782*GR;
    hs_S=2.9-782*GR_S;
    hs_AMSR= 2.9-782*GR_AMSR;
    
    figure
    grid on
    hold on;
    
    plot(delta_density,hs-hs_AMSR,'r');
    plot(delta_density,hs_S-hs_AMSR,'g');


%    plot(hs_AMSR*ones(1,30),delta_di,'k');
    title('Snow depth retrieval difference vs snow density');
       xlabel('Snow density [kg/m^3]');
       ylabel('Snow depth retrieval difference [cm]');
    legend('without','with');
    hold off 

    
    end
 

    if n==3
      
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
 
    GR = (TB_FY_37G(1,:)-TB_FY_19G(1,:))./(TB_FY_37G(1,:)+TB_FY_19G(1,:));    
    GR_S = (TB_FY_37G_S(1,:)-TB_FY_19G_S(1,:))./(TB_FY_37G_S(1,:)+ ...
                                                TB_FY_19G_S(1,:));
GR_AMSR = (output(1,9)-output(1,5))/(output(1,9)+output(1,5));


    hs= 2.9-782*GR;
    hs_S=2.9-782*GR_S;
    hs_AMSR= 2.9-782*GR_AMSR;
    
    figure
    grid on
    hold on;
    
   plot(delta_pci,hs-hs_AMSR,'r');
   plot(delta_pci,hs_S-hs_AMSR,'g');
      title('Snow depth retrieval difference vs snow density');
    xlabel('Snow correaltion length [mm]');
    ylabel('Snow depth retrieval difference [cm]');
    legend('without','with');
    hold off 

    
    end
 

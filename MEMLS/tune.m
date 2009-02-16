
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             K               kg/m^-3 cm 
%%      num,   T          W   ,  roi   ,    di ,     pci        sal  type
%%  par = [1 , Tsea ,   0.0026, roi_ice,   di_ice,   0.35,      7.5,  1;
%%         2 , Ts_snow ,0.0000, roi_snow2, di_snow2, pci_snow2, 0,    0;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = tune (n)
% Tb under different profiles on same plot  
% n: 1/2/3 snow density/correlation length/depth

  
par_in = load('input_FY') ;
delta_density = 30;
delta_pci = 0.01;
delta_di = 5;

output = zeros(1,10);
TB_MY= zeros(1,10);
TB_FY=zeros(1,10);

freq = [6.9,10.7,18.7,23.8,36.5];

cd Data;
output = average_T('area.58s.sort',2);
cd ..;


%% ################### PLOT ############################### 
figure
grid on
hold on

plot(freq,output(1:2:10),'--ms', 'LineWidth',1.5);%V channel
plot(freq,output(2:2:10),'--gs', 'LineWidth',1.5);%H channe


niter = 5;
rval = linspace(0,1,niter);
gval = linspace(0,1,niter);

if n ==1 
for i = 1:niter
  par_in(3,4) = par_in(3,4) + delta_density;
  
  par_out =  homo(par_in,10,2);
  par_out = par_out';
  TB_FY=fw_fun2(par_in,2,0.2);
  
  %      set(gcf,'Color',[1,0.4,0.6])
  plot(freq,TB_FY(1:2:10),'Color',[rval(i),0.8,0.8]);
  plot(freq,TB_FY(2:2:10),'Color',[0.8,gval(i),0.8]);
  
  xlabel('Frequency [GHz]');
ylabel('Brightness Temperature [K]');
title('FY: Tb vs. Frequency in Jan of area 58 at various snow density');
legend('Tbv','Tbh','Tbv_a','Tbh_a',2);

  
  fid = fopen('output_FY', 'wt');
  fprintf(fid,'%1d %3.2f %3.5f %4.4f %4.4f  %4.4f  %4.4f %1d\n',par_out');
  fclose(fid);
end
hold off
end


if n ==2 
for i = 1:niter
  par_in(3,6) = par_in(3,6) + delta_pci;
  
  par_out =  homo(par_in,10,2);
  par_out = par_out';
  TB_FY=fw_fun2(par_in,2,0.2);
  
  %      set(gcf,'Color',[1,0.4,0.6])
  plot(freq,TB_FY(1:2:10),'Color',[rval(i),0.8,0.8]);
  plot(freq,TB_FY(2:2:10),'Color',[0.8,gval(i),0.8]);
  
  fid = fopen('output_FY', 'wt');
  fprintf(fid,'%1d %3.2f %3.5f %4.4f %4.4f  %4.4f  %4.4f %1d\n',par_out');
  fclose(fid);
  xlabel('Frequency [GHz]');
ylabel('Brightness Temperature [K]');
title('FY: Tb vs. Frequency in Jan of area 58 at various snow correlation length');
legend('Tbv','Tbh','Tbv_a','Tbh_a',2);

end
hold off
end

if n ==3
for i = 1:niter
  par_in(3,5) = par_in(3,5) + delta_density;
  
  par_out =  homo(par_in,10,2);
  par_out = par_out';
  TB_FY=fw_fun2(par_in,2,0.2);
  
  %      set(gcf,'Color',[1,0.4,0.6])
  plot(freq,TB_FY(1:2:10),'Color',[rval(i),0.8,0.8]);
  plot(freq,TB_FY(2:2:10),'Color',[0.8,gval(i),0.8]);
  
  fid = fopen('output_FY', 'wt');
  fprintf(fid,'%1d %3.2f %3.5f %4.4f %4.4f  %4.4f  %4.4f %1d\n',par_out');
  fclose(fid);
  xlabel('Frequency [GHz]');
ylabel('Brightness Temperature [K]');
title('FY: Tb vs. Frequency in Jan of area 58 at various snow depth');
legend('Tbv','Tbh','Tbv_a','Tbh_a',2);

end
hold off
end


% %%%%%%%%%%%%CONSTANT%%%%%%%%%%%%%%%%%%%%%%%%%
% Tsea=271.35;      %[K] 
% ks=0.31;
% k0=2.034;         %[W/m/K]
% beta=0.117        %[W/m/per mil] 

% %%%%%%%%%%%Variable to be set%%%%%%%%%%%%%%%
% Ts_snow = ; 
% roi_snow = ;
% di_ice = ;
% di_snow1 = 
% di_snow2 =
% pci_snow = ;
% sal_ice = ;
% sal_2 = ;
% Range of validity of surface tempurature of snow surface

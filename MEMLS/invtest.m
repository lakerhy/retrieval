function[p,p2,p3]= invtest()
input = load('input_FY');
    V = 2;
Ts_MY = input(2,1);
di_snow = input(2,3);
roi_snow = input(2,2);
pci_snow = input(2,4);
di_ice = input(1,3);
sal = input(1,5);
    
    %emissivity_MY=Ts_MY./Tb_MY 
p = [Ts_MY; V;di_snow;roi_snow;pci_snow;di_ice;sal ];    
TB_MY=fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1),2);
     
    figure
    grid on
    hold on
    plot(1:5,TB_MY(1:2:10),'-rs');
    plot(1:5,TB_MY(2:2:10),'-gs');

p2 = inv_funtest(TB_MY,2);
%TB_MY2=fw_fun2(p2);
TB_MY2=fw_fun2(p2(1,1),p2(2,1),p2(3,1),p2(4,1),p2(5,1),p2(6,1),p2(7,1),2);
     
     hold on
    plot(1:5,TB_MY2(1:2:10),'--ms');%,'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
    plot(1:5,TB_MY2(2:2:10),'--bs');

legend('Tb_v','Tb_h','Tb_v_i','Tb_h_i'); 
hold off;
p2
p
p3= p2(:,3)-p

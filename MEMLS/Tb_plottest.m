[p1] = retrieval(1)
p=p1';


    TB_MY=fw_fun2(p(1,1),p(2,1),p(3,1),p(4,1),p(5,1),p(6,1),p(7,1));
    
    
    
 cd Data;
    output = average_T('area.33s.sort');
    cd ..;
    
        figure
    grid on
    hold on
    plot(1:5,TB_MY(1:2:10),'-rs');%,'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
    plot(1:5,TB_MY(2:2:10),'--gs');
    %V channel
    plot(1:5,output(1:2:10),':ms');
    %H channel
    plot(1:5,output(2:2:10),'-.bs');
    title('MY: Brightness temperature vs. Frequency in Jan of location.33 ');

    legend('Tbv','Tbh','Tbv_a','Tbh_a',2);
    hold off
inv58 = load('58dayinv');
depth = inv58(:,2);
jan = [1:29];
feb = [30:57];
mar = [58:88];
apr = [89:118];
may = [119:149];

figure
hold on 
grid
plot(jan,depth(jan),'r+');
plot(feb,depth(feb),'g+');
plot(mar,depth(mar),'b+');
plot(apr,depth(apr),'k+');
plot(may,depth(may),'m+');
hold off

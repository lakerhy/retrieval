function[] = loc (filename,range)

cd ../MEMLS/Data
buffer = load(filename);
f1V = buffer(:,6); % 6
f1H = buffer(:,7);% 17
f2V = buffer(:,8);
f2H = buffer(:,9);
f3V = buffer(:,10); % 6
f3H = buffer(:,11);% 17
f4V = buffer(:,12);
f4H = buffer(:,13);
f5V = buffer(:,14); % 6
f5H = buffer(:,15);% 17

lat = buffer(:,16); % 16
lon = buffer(:,17);% 17

cd ../../stat

num = size(f1V);
f1pol = f1V-f1H;
f2pol = f2V-f2H;
f3pol = f3V-f3H;
f4pol = f4V-f4H;
f5pol = f5V-f5H;

figure 
subplot(2,3,1)
plot(lon,lat,'wd','MarkerFaceColor','r','MarkerSize',6);
xlabel('Longitude');
xlabel('Latitude');
title(['location distribution of ',filename]);

lon1 = range(1);
lon2 = range(2);
lat1 = range(3);
lat2 = range(4);

lat_sub = find(lat<=lat2&lat>lat1);
lon_sub = find(lon<=lon2&lon>lon1);
sub = intersect(lat_sub,lon_sub);

subplot(2,3,2);
plot(1:size(sub),f1pol(sub),'wd','MarkerFaceColor','r','MarkerSize',6);
ylabel('Pol.')
legend('6.9G');

subplot(2,3,3);
plot(1:size(sub),f2pol(sub),'wd','MarkerFaceColor','r','MarkerSize',6);
ylabel('Pol');
legend('11G');

subplot(2,3,4);
plot(1:size(sub),f3pol(sub),'wd','MarkerFaceColor','r','MarkerSize',6);
ylabel('Pol');

legend('19G');

subplot(2,3,5);
plot(1:size(sub),f4pol(sub),'wd','MarkerFaceColor','r','MarkerSize',6);
ylabel('Pol');
legend('24G');

subplot(2,3,6);
plot(1:size(sub),f5pol(sub),'wd','MarkerFaceColor','r','MarkerSize',6);
ylabel('Pol.');
%xlabel('sample number ');
legend('37G');


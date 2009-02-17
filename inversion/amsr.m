function[Tb_amsr] = amsr (filename)
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

f1Vmean = mean(f1V);
f1Hmean = mean(f1H);
f2Vmean = mean(f2V);
f2Hmean = mean(f2H);
f3Vmean = mean(f3V);
f3Hmean = mean(f3H);
f4Vmean = mean(f4V);
f4Hmean = mean(f4H);
f5Vmean = mean(f5V);
f5Hmean = mean(f5H);

Tb_amsr = [f1Vmean,f1Hmean,f2Vmean,f2Hmean,f3Vmean,f3Hmean,f4Vmean,f4Hmean,f5Vmean,f5Hmean];
cd ../../inversion

function ta=transang(ia,e1,e2)
%conputes the transmission angle using Snellius
ta=asin((sqrt(real(e1)) ./sqrt(real(e2))) .*sin(ia));

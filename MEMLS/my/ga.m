function[g] = ga (freq,epsi)

   c = 2.99793;
lamd=c/(10*freq);
g= 2*pi/lamd*sqrt(epsi);

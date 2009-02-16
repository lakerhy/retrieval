% beregner dielektrisitetskonstanten af våd sne

function e_ws=e_wet_snow(k,Dens,wc)

   fws=k./20.94;
   denws=Dens./1000.0;
   wc=wc.*100.0;
   A1ws=1.0;
   A2ws=1.0;
   B1ws=0.0;
   Aws=1.0+1.83.*denws+0.02.*A1ws.*(wc.^1.015)+B1ws;
   Bws=0.073.*A1ws;
   Cws=0.073.*A2ws;
   f0=9.07; %relaxation frequency

   ews1=Aws+((Bws.*wc)./(1.0+((fws./f0).^2)));
   ews2=(Cws.*(fws./f0).*wc)./(1.0+((fws./f0).^2));
   e_ws=ews1+ews2*i;

end

function [rho,dzo]=snow_melt(rho,dzi,Qtot,Mr,Cwf)

%;if (dzi<=0.001 | rho>=920) 
%; dzo=0.001; rho=920; 
%;return end

Lf=0.334e6;
dzo=dzi;

%tykkelse
if (Qtot<0) dzo=dzi; end
if (dzi>Qtot./(Lf.*rho) & Qtot>0.001) dzo=dzi-Qtot./(Lf.*rho); end
if (dzi<=Qtot./(Lf.*rho) & Qtot>0.001) dzo=0.001; end

%massefylde
f=(dzi-dzo)./dzi;
drhow=(Cwf+Mr)./dzi;
rho=(rho+drhow)./(1-f);


end

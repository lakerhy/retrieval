function [wc,dz]=layer_melt(Cw,rho,d)
Cw_max=rho.*d;
if (Cw>Cw_max) wc=1.0; dz=0; break; end
wc=(Cw)./Cw_max;
dz=(Cw_max-Cw).*d./Cw_max;
end

function U=water_flux(wc,time,gs,dens)
%brun et al. 1989
Sw=wc;
%k=8e-10; %snow permeability from Albert, 1996
k=0.077.*gs^2.*exp(-7.8.*dens./1000.); %Gallee&duynkerke,1997,JGR102(D12),13813-13824.
g=9.8; %gravitational acceleration
Swi=0.04; %irreducible water saturation, 4% Jordan et al 1999, 3% Albert&Krajeski,1998
rhow=1000.; %density of water
myw=1.79e-3; %viscosity of water

U=0;
if (Sw>Swi)
  U=time.*(rhow.*g.*k./myw).*((Sw-Swi)./(1-Swi)).^3;
end
if (U>Sw-Swi) U=Sw-Swi; end

end

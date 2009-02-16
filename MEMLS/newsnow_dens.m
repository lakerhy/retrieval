% after Jordan et al. 1999 JGR 104 C4 7785-7806, eq. 14
% beregn nyfalden snes massefylde som funktion af luft temperatur og vindhastighed
function rho=newsnow_dens(Tair,u)
 if (Tair>260.15)
   rho=500.0*(1.0-0.951.*exp(-1.4.*((278.15-Tair).^-1.15)-0.008.*u.^1.7)); 
 else
   rho=500.0.*(1.0-0.904.*exp(-0.008.*u.^1.7));
 end %endif
end %endfunction

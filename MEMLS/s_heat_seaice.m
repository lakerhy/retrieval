function csi=s_heat_seaice(rho,S,T)
%computes the specific heat of sea ice given the SH of pure ice Makshtas (1998)
cpi=2000.0; %is actually a function of temp
csi=cpi+(1.715e7 .*S) ./(rho .*(T-273.16).^2);

end

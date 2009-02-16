function pec=d2pec_matz(d,dens)
%Matzler, 2002, J. Glaciol.
%d[m]
%dens[kg/m3]
%pec[mm]
v=dens ./927.0;
pec=500.0 .*d .*(1.0-v);
end

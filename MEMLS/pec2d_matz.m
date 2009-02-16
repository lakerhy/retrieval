function d=pec2d_matz(pec,dens)
%Matzler, 2002, J.Glaciol.
v=dens ./927.;
d=(pec .*2.0) ./(1000.0 .*(1.0-v));
end

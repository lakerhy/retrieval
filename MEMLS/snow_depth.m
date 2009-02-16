function d=snow_depth(d,density,T,load,t)
%jordan, 1991
mu0=3.6e6; %N/m2
c3=1;
c4=1;
c5=0.08; %1/K
c6=0.021; %m3/kg

if (density>150.0)
   c3=exp(-0.046.*(density-150.0));
end %endif

   meta=2.778e-6 .* c3 .* c4 .* exp(-0.04.*(273.15-T));

   overb=(load./mu0).* exp(-c5.*(273.15-T)).*exp(-c6.*density);

   CR=meta+overb;

d=d-d.*CR.*t;

end

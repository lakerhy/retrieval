function F=shortwave_energy1(type,sw_i,dz,n,rho,pec,time)

%computes short wave energy input to each layer
%albedo and extinction depend on the type

%albedo=shortwave_albedo(type(1));
[ac,albedo]=crocus_albedo(pec2d(pec(1)),rho(1));
if (type(1)>2 || rho(1)>700.) albedo=0.72; end %endif

F(1:n)=0;

energy_a=(sw_i.*(1-albedo));

extinction=shortwave_extinction(type,n,rho,pec);

for i=1:n
 energy_b=energy_a.*exp(1).^(-extinction(i).*dz(i));
 F(i)=energy_a-energy_b;
 energy_a=energy_b;
 if (energy_a<1e-5) break; end
end

F=time.*F;

end

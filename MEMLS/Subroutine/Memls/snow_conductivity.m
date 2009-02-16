function ks=snow_conductivity(rho,Ts)

%computes the thermal snow conductivity Makshtas (1998) eq. 14

ks=2.845e-6 .* rho.^2 + 2.7e-4 .* 2.^((Ts-233)./5);
end

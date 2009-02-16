function ki=por_ice_conductivity(S,T,rho)

%computes the thermal conductivity of porous (sea) ice Makshtas (1998) eq. 9-10
ka=0.024;
va=1.0-(rho./920.0);

if S==0.0
 kbi=2.22.*(1.0-0.0159.*(T-273.16));
else
 kbi=sal_ice_conductivity(S,T);
end

f=ka./kbi;

ki=kbi.*(1.0+0.5.*f-va.*(1.0-f))./(1.0+0.5.*f+0.5.*va.*(1.0-f));

end

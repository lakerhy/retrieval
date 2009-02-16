function ki=sal_ice_conductivity(S,T)

%computes the conductivity of saline ice after Makshtas (1998) eq. 12
if T<272 
  ki=2.03+(0.12.*S)./(T-273.0);
%  ki=9.828.*exp(-0.0057.*(273.0+(T-273.15)))+((0.13.*S)./(T-273.15));
else
  ki=2.0;
end %endif

end

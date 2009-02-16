function sh=sensible_heat(u,Ta,T0,time)

rhoa=1.25; %density of air see eg Garratt (1992) p. 284
cpa=1005.0; %specific heat of air see eg Garratt (1992) p. 284
Cs=1.2e-3; %transfer coef. see Maykut 1986 p. 415

sh=time.*rhoa.*cpa.*Cs.*u.*(Ta-T0);

end

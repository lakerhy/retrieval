function lh=latent_heat(u,r,Ta,T0,P,time)

rhoa=1.25; %dens of air see garratt 1992 p. 284
L=2.502e6; %latent heat of vaporisation see garratt 1992 p. 284
Ce=0.55e-3; %bulk transfer coef see Maykut 1986 p. 415
%es(a or 0) are the partial water vapour pressure in the air and at the surface respec.
%see garratt 1992 p. 284
esa=6.112.*exp(1).^(17.67.*(Ta-273.15)./(Ta-29.65));
es0=6.112.*exp(1).^(17.67.*(T0-273.15)./(T0-29.65));
%r is the relative humidity at measurement height
r=r./100;
%P is the air pressure

lh=time.*0.622.*rhoa.*L.*Ce.*u.*((r.*esa)-es0)./P;
end

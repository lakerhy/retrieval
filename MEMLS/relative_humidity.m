function rh=relative_humidity(Ta,Td)

Ta=Ta-273.15;
Td=Td-273.15;

vd1=5.018+0.32321.*Ta+8.1847e-3.*(Ta.^2)+3.1243e-4.*(Ta.^3);
vd2=5.018+0.32321.*Td+8.1847e-3.*(Td.^2)+3.1243e-4.*(Td.^3);
rh=abs(100.0.*vd2./vd1);

end

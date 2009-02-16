function c=s_heat_salmix(s,T,M)

if (M==0) 
    return 
end
%computes specific heat of saline snow or ice mixture
Lw=0.334e6; %latent heat of fusion
cpi=2113.0; %spec. heat of pure ice
%cbr=4217.0; %spec. heat of water (brine)
cbr=s_heat_sw(Sb(T-273.15),T-273.15,0); %sp.heat of saline water
[dummy,da]=swstate(Sb(T-273.15),T-273.15,0);
db=1000+da; %density of brine
vb=Vb(T-273.15,s); %volume of brine
dvb_dT=Vb(T-273.15,s)-Vb(T-274.15,s);
%s=0.001.*s;
%sb=0.001.*Sb(T-273.15); %salinity of brine
Mi=M-vb.*db; %mass of pure ice
Mb=vb.*db; %mass of brine

c=cpi.*(Mi./M)+cbr.*(Mb./M)+Lw.*Mb.*dvb_dT;

end

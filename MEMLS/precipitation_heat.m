function Qr=precipitation_heat(Mr,Ta,T0,time)
%energiinput fra regn per tidsenhed
%brun et al. 1989
cpw=4217.0; %specific heat of water
Qr=Mr.*cpw.*(Ta-T0);
end

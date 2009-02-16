function pec=depth_hoar(pec,T,dTdz,rho,time)

dTdz=abs(dTdz);

%udviklingen af depth hoar i t?r sne
%marbouty, 1980
if (rho >= 400) 
    return; 
end

a0=pec./400.0; %correlation length til optisk diameter

f=0.0;
h=0.0;
g=0.0;

%function f
if (T>266) f=0.7+0.043.*(273-T); end
if (T<=266 & T>251) f=1-0.053.*(266-T); end
if (T<=251 & T>233) f=0.2-0.011.*(251-T); end
if (T<=233) f=0; end

%function h
if (rho<150) h=1; end
if (rho>=150 & rho<400) h=1-0.004.*(rho-150); end
if (rho>=400) h=0; end

%function g
if (dTdz<0.15) g=0.0; end
if (dTdz>=0.15 & dTdz<0.7) g=1.82.*(dTdz-0.15); end
if (dTdz>=0.7) g=1.0; end

o=1.04e-9.*time;

a=a0+(f.*g.*h.*o);

pec=400.*a;

end

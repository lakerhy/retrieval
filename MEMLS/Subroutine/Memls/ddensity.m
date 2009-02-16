function rho=ddensity(rho,load,time)
%Sturm & Holmgren, 1998 eq. 3.
k=0.02;
ny0=8.5e6;
%viskositet
rho=rho+((load.*rho.*time)./(ny0.*exp(k.*rho)));

end

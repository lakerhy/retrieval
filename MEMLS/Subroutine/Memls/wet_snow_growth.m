function pec=wet_snow_growth(pec,wc,time)

%brun 1989, ann. glac. 13,22-26
%wc=100.*wc;
v00=1.28e-8;
v11=4.22e-10;
r=pec./0.8;
v0=1.33.*3.14.*(r.^3);
vt=time.*(v00+v11.*(wc.^3));
v=v0+vt;
r=((3.*v)./(4.*3.14)).^0.333;
pec=r.*0.8;

end

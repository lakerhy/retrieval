function grad=dTdz(T1,T2,dz)
if (dz<0.001) grad=(T1-T2)./0.001; disp('thin top layer: dTdz\n');
else  grad=(T1-T2)./dz; end
end

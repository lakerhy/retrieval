function ss=iborn_s2p_bc(e1,e2,eeff,v,k,pcc)
% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
%backscattering coefficient of a collection of spherical
%particles with correlation length pcc using improved born approximation
R=3. .*pcc ./4.;
ss=(R.^3 .*k.^4 ./3) .*v .*(1-v) .*abs(((e2-e1) ...
.*(2 .*eeff+e1)) ./(2 .*eeff+e2)).^2;

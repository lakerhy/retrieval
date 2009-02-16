function ac=iborn_s2p_ab(e1,e2,eeff,v,k)
% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
%absorption coefficient of a collection of spherical
%particles with correlation length pc using improved born approximation
ac_s= v.*k.*abs(imag(e2)).*abs((2.0 .*eeff + e1) ./(2.0 .*eeff + e2)).^2;
ac_b= 2.0 .*(1-v).*k.*abs(imag(sqrt(e1)));
ac=ac_s+ac_b;

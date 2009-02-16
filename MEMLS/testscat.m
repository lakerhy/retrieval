function [ss,sphe]=testscat(e1,e2,eeff,v,k,pcc)
% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
%scattering coefficient of a collection of spherical
%particles with correlation length pc using improved born approximation
ss=(3 .*pcc.^3 .*k.^4 ./32) .*v .*(1-v) .*abs(((e2-e1) ...
.*(2 .*eeff+e1)) ./(2 .*eeff+e2)).^2;
pci=pcc.*1000.0;
epseff=eeff;
vfi=v;
eice=e2;
sphe = (3/32).*(0.001.*pci).^3.*k.^4.*vfi.*(1-vfi).*abs((2.*epseff+1).*(eice-1)./(2.*epseff+eice)).^2;

endfunction

%function [eeff,epseff]=testscat(e1,e2,v)

% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
% Polder/VanSanten mixing formulae for spheical inclusions
%effective dielectric constant of medium consisting of e1 and e2
%e1: dielectric constant of background
%e2: dielectric constant of sherical inclusion
%v: fraction of inclusions

%eeff=0.25 .*(2 .*e1-e2+3 .*v .*(e2-e1)+sqrt((2 .*e1-e2+3 ...
%.*v .*(e2-e1)).^2 +8 .*e1 .*e2));

%eice=e2;
%vfi=v;

%epseff = (2-eice+3.*vfi.*(eice-1)+ sqrt((2-1+3.*vfi.*(eice-1)).^2+8.*eice))./4;
%epseff = (2-eice+3.*vfi.*(eice-1)+ sqrt((2-eice+3.*vfi.*(eice-1)).^2+8.*eice))./4;

%endfunction

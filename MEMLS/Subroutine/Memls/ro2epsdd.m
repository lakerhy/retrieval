 function [epsi,epsii] = ro2epsdd(roi,Ti,freq,graintype)

 % calculates the dielectric permittivity from
 % density for dry snow.
 %
 % [epsi,epsii] = ro2epsd(roi,Ti,freq)
 % epsi: real part of dielectric permittivity
 % epsii: imaginary part of dielectric permittivity
 % roi: density
 % Ti: snow temperature in Kelvin
 % freq: frequency
 %
 % Version history:
 % 1.0 wi 15.7.95
 % 2.0 wi 12.11.97 enhanced with Polder and van Santen Equations (see Polder.m)
 %
 % Uses:
 % epsice, epsr, polder
 %
 %
 %
 % Copyright (c) 1997 by the Institute of Applied Physics,
 % University of Bern, Switzerland

 eice = epsice(Ti,freq);

 epsi = epsr(roi);

 % imaginary part after Tiuri 84
 %epsii = eice.*(0.52.*roi + 0.62.*(roi.^2));

 if graintype==1 %for the case of snow with the empirical formula
 % imaginary part after Polder and van Santen 1946 (Effective-Medium Approx)
 ei=3.185;
 eice = epsice(Ti,freq);
 eps2=ei;
 f = roi ./ 0.917;
 N = length(roi);
 A = zeros(N,1);
 A = A + 1/3;
 for i=1:N
 if f(i) < 0.55
 A(i) = 0.476 - 0.64 * f(i);
 end
 if f(i) <= 0.333
 A(i) = 0.1 + 0.5 * f(i);
 end
 end
 epseff=epsr(roi);
 A3 = 1 - 2 .* A;
 ea = epseff .* (1-A) + A;
 ea3 = epseff .* (1-A3) + A3;
 K1 = (ea ./ (ea+A .* (eps2-1))).^2;
 K3 = (ea3 ./ (ea3+A3 .* (eps2-1))).^2;
 Ksq = (2 .* K1 + K3) ./ 3;
 epsii = sqrt(epseff) .* eice .* Ksq .* f;

 elseif graintype==2 %case of small spherical scatterers
 ei=3.185;
 eice = epsice(Ti,freq);
 f = roi ./ 0.917;
 eps1=1;
 eps2=ei;

 epseff=(2*eps1-eps2+3*(f)*(eps2-eps1)+sqrt((2*eps1-eps2+3*(f)*(eps2-eps1)).^2+8*eps1*eps2))/4;
 Ksq = ((2*epseff+eps1)./(2*epseff+eps2)).^2;
 epsii = sqrt(epseff) .* eice .* Ksq .* f;

 elseif graintype==3 %case of thin shells
 ei=3.185;
 eice = epsice(Ti,freq);
 f = roi ./ 0.917;
 epseff=eps1+(f*(eps2-eps1)*(2+eps1/eps2))/(3-f*(1-eps1/eps2));
 eps1=1;
 eps2=ei;
 Ksq = (2/3) + (eps1.^2)./(3*(eps2.^2));
 epsii = sqrt(epseff) .* eice .* Ksq .* f;
 else
 disp('Problem with the assumption of the mixing formula')
 end

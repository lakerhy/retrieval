function [epsi, epsii] = te(sal,epsi,epsii)

% test function for sea ice dielectrics

  aepsi=3.15;
  aepsii=0.2;

  epsi=epsi-epsi.*ceil(sal)+aepsi.*ceil(sal);
  epsii=epsii-epsii.*ceil(sal)+aepsii.*ceil(sal);

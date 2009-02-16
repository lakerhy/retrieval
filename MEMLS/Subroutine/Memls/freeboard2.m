function [sfb,dlay]=freeboard2(d,Dens)

%beregner snefrihøjden over vandoverfladen og forsinkelsen i tid
%for en radar 1m over vandoverfladen.

  dlay=0;
  a=1000-Dens(2);
  b=d(2);
  c=d(1).*Dens(1);
  d=Dens(2);
  mass=sum(d.*Dens);
  height=sum(d);
  sfb=((a.*b)-c)./(d+a);
  dlay=(1.0-sfb)./3e8;

endfunction

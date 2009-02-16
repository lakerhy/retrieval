function [sfb,dlay]=freeboard(d,Dens)

%beregner snefrih?jden over vandoverfladen og forsinkelsen i tid
%for en radar 1m over vandoverfladen.

  dlay=0;
  mass=sum(d.*Dens);
  height=sum(d);
  sfb=height.*(1.0-(mass./(height.*1000.0)));
  dlay=(1.0-sfb)./3e8;

end %endfunction

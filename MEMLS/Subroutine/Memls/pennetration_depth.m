%computes the penetration depth, using radiative transfer
%n: number of layers
%vol: volume scattering coefficient
%surf: surface scattering coefficient
%loss: loss coefficient
%trans: transmissivity of layer surface

function [dpt,x]=pennetration_depth(n,d,loss,trans)

inten=1.0;
einv=0.367879;
penne=0.0;
dpt=0.0;
x=0.0;

mtrans=trans;
mtrans(1)=1.0;

inten(1:n)=1.0;

for idx=1:n
  inten(idx)=prod(1./(loss(1:idx))) .*prod(mtrans(1:idx));
  penne(idx)=sum(d(1:idx));
  dpt=penne(idx);
  if inten(idx) <= einv
    if idx==1
      slope=(1.0-inten(idx))./d(idx);
      x=(1.0-einv)./slope;
    else
      slope=(inten(idx-1)-inten(idx))./d(idx);
      x=penne(idx-1)+(inten(idx-1)-einv)./slope;
    endif
    break;
  endif
endfor

endfunction

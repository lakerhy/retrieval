function loss=extloss(d,ang,ke)
%computes the loss coefficient in a layer with extinction, ke, and thickness, d.

loss=exp(ke .*d ./cos(ang));

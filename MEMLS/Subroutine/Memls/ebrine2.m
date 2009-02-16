function ebr=ebrine2(freq,T,Si)

% permittivity and loss of brine Ulaby et al. 1986 E64a
% freq: em frequency

volb=Vb(T,Si);
salb=Sb(T);
N=Nsw(salb);
conb=condbrine(T,N);
relax=relaxt(T,N);
eb0=epsib0(T,N);

f=freq .*1.0e9;
epsiwoo=4.9;
e0=8.854e-12;

eb=epsiwoo+((eb0-epsiwoo) ./(1.+((relax .*f) .^2)));
ebi=((relax .*f .*(eb0-epsiwoo)) ./(1.+((relax .*f) .^2)))+(conb ./(6.28 .*f .*e0));
ebr=eb+i*ebi;

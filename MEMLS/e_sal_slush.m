%beregner dielektricitetskonstanten af saltholdigt sjap
function ess=e_sal_slush(Dens,wc,S,TK,freq)

%Dens: massefylde i kg/m3
%wc: vandindhold (volumen)
%S: saliniteten i ppt
%TK: temperaturen i kelvin
%freq: frekvensen i GHz

eair=1.0+0.0i;
T=TK-273.15;

%beregn saliniteten af væsken, salt findes kun i væsken
Sv=(Dens .*S) ./(1000 .*wc);

%beregn den dielektriske konst. af luft og væske (baggrund)
N=Nsw(Sv);
conb=condbrine(T,N);
relax=relaxt(T,N);
eb0=epsib0(T,N);
f=freq .*1.0e9;
epsiwoo=4.9;
e0=8.854e-12;
eb=epsiwoo+((eb0-epsiwoo) ./(1.+((relax .*f) .^2)));
ebi=((relax .*f .*(eb0-epsiwoo)) ./(1.+((relax .*f) .^2)))+(conb ./(6.28 .*f .*e0));
ebr=eb+i*ebi;
eba=eice_s2p(eair,ebr,wc);

%beregn den dielektriske konstant af ren is
epui=e_ice(TK,freq);

%beregn den dielektriske konstant af slush
v=Dens ./(1000 .*wc + 920 .*(1-wc));
ess=eice_s2p(eba,epui,v);

endfunction

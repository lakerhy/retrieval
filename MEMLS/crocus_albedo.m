function [ac,agk]=crocus_albedo(gs,dens)

%computes the albedo using CROCUS without zenit angle correction and Greuell & Konzelmann, 1994.
alpha1=1-1.58.*sqrt(gs);
alpha2=1-15.4.*sqrt(gs);
alpha3=346.3.*gs-32.31.*sqrt(gs)+0.88;
ac=0.606.*alpha1+0.301.*alpha2+0.093.*alpha3;

agk=0.58+(dens-920).*-4.3548e-4;
end

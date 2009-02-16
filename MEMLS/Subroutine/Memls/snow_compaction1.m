function d=snow_compaction1(d,load,dtime,density,Tn,Tmelt)

%computes the compaction after brun et al. 1989 perhaps corrected for errors!
density=density./1000.0;
dtime=dtime./3600.0;
fd=0.4;

func=((6e4) .*e.^((23.0 .*density)-(0.1 .*(Tn-Tmelt))))./(1.0-fd);
%compaction of unit layer depth
dd=(load.*dtime)./func;
if ((1-dd) < density./0.92)
% d=(density./0.92).*d;
 d=0.001;
 density=0.92;
else
 %compute new density
 density=(1+dd).*density;
 d=(1-dd).*d;
end

end

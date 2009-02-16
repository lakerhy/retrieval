function Ti=dry_temperature(Tabs,totalq,varmef,massef,tyk)

%beregner lagets temperatur udfra energibalancen, den maksimale ?ndring over 10min er 0.7C

if (varmef.*massef.*tyk == 0) 
 Ti=Tabs;
 disp('varmef | massef | tyk =0');
 return;
end %endif

dT=(totalq./(varmef.*massef.*tyk));
if (totalq./(varmef.*massef.*tyk)) > 0.5 dT=0.5; end %endif
if (totalq./(varmef.*massef.*tyk)) < -0.5 dT=-0.5; end %endif

Ti=Tabs+dT;

if (Ti>271.0) Ti=271.0; end %endif
end %endfunction

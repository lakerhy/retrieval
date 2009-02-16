function [snow,ice]=thickness(n,d,typ)

ice=0;
snow=0;

for i=1:n

if (typ(i)<=2) snow=snow+d(i); end %endif
if (typ(i)>2) ice=ice+d(i); end %endif

end %endfor

end %endfunction

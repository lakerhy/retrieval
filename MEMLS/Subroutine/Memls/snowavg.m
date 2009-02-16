function [adens,apcc,asT,asal,aiT,isT,pcc1]=snowavg(n,typ,d,dens,pcc,T,sal)
%average properties of the snow
snow=(typ<=2);
if (sum(snow)<1) 
    adens=0;
    apcc=0;
    asT=0;
end %endif

%average properties of the ice
ice=(typ>2);
ice
if (sum(ice)<1)
    asal=0;
    aiT=0;
    isT=0;    
end %endif

if (sum(snow)>=1)
adens=sum(snow.*d.*dens)./sum(snow.*d);
apcc=sum(snow.*d.*pcc)./sum(snow.*d);
asT=sum(snow.*d.*T)./sum(snow.*d);
end %endif

o=1;
%if (sum(ice)>=1)
asal=sum(ice.*d.*sal)./sum(ice.*d);
aiT=sum(ice.*d.*T)./sum(ice.*d);
for i=1:n
    if (ice(i)>0)
        o=i;
        break;
    end
end

isT=T(o);
pcc1=pcc(o);

%end %endif
end %endfunction

function [n,num,T,type,Dens,d,pcc,sal,wc,rms,si]=del_layer(n,num,T,type,Dens,d,pcc,sal,wc,rms,si)

%slet lag hvis det bliver for tyndt <1mm.

%l?s parametre i temp variabel
tnum=num; tT=T; ttype=type; tDens=Dens; td=d; tpcc=pcc; tsal=sal; twc=wc; trms=rms; tsi=si;

idx=0;

for i=1:n
 if (td(i)>0.001)
  idx=idx+1;
  num(idx)=idx; T(idx)=tT(i); type(idx)=ttype(i); d(idx)=td(i); pcc(idx)=tpcc(i);
  sal(idx)=tsal(i); wc(idx)=twc(i); rms(idx)=trms(i); si(idx)=tsi(i); Dens(idx)=tDens(i);
 end %endif
end %endfor
num=num(1:idx); T=T(1:idx); type=type(1:idx); d=d(1:idx); pcc=pcc(1:idx); 
sal=sal(1:idx); wc=wc(1:idx); Dens=Dens(1:idx); rms=rms(1:idx); si=si(1:idx);

% n=max(size(d));
n=length(d);

end %endfunction

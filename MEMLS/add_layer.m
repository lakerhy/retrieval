function [n,num,T,type,Dens,d,pcc,sal,wc,rms,si]=add_layer(n,num,T,type,Dens,d,pcc,sal,wc,rms,si,Ms,Ta,u)

%nye snelag har massefylden: 
nDens=newsnow_dens(Ta,u);
dns=Ms./nDens;
if (dns>0)
 %l?s parametre i temp variabel
 tnum=num; tT=T; ttype=type; tDens=Dens; td=d; tpcc=pcc; tsal=sal; twc=wc; trms=rms; tsi=si;
 %s?t f?rste lag
 num(1)=1; T(1)=Ta; type(1)=1; Dens(1)=nDens; d(1)=dns; pcc(1)=0.06; 
 sal(1)=0.0; wc=0.0; rms=0.015; si=0.0;
 for i=2:n+1
   num(i)=i; T(i)=tT(i-1); type(i)=ttype(i-1); Dens(i)=tDens(i-1); d(i)=td(i-1);
   pcc(i)=tpcc(i-1); sal(i)=tsal(i-1); wc(i)=twc(i-1); rms(i)=trms(i-1); si(i)=tsi(i-1);
 end %endfor
end %endif

% n=max(size(num));
n=length(num);

%hvis bundlaget er blevet over 6cm tykt dannes et nyt lag.
%havis gror fra bunden:
if (d(n)>0.06)
  d(n+1)=d(n)-0.05;
  d(n)=0.05;
  num(n+1)=n+1; T(n+1)=0.5.*(T(n)+271.35); type(n+1)=type(n); Dens(n+1)=Dens(n); pcc(n+1)=pcc(n);
  sal(n+1)=sal(n); wc(n+1)=0; rms(n+1)=0; si(n+1)=si(n);
  n=n+1;
end %endif

end %endfunction

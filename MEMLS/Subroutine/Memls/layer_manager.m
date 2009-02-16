function [n,num,T,type,Dens,d,pcc,sal,wc]=layer_manager(n,num,T,type,Dens,d,pcc,sal,wc,Ms,Ta)

%idx3=(d>0.001);
%T=T(idx3); type=type(idx3); Dens=Dens(idx3); d=d(idx3);
%pcc=pcc(idx3); sal=sal(idx3); wc=wc(idx3);
%n=max(size(T));
%for i=1:n
% num(i)=i;
%end

%læs parametre i temp variabel
tnum=num; tT=T; ttype=type; tDens=Dens; td=d; tpcc=pcc; tsal=sal; twc=wc;

%num
%T
%Dens
%d

%idx3=(d>0.001);
%yy=d(idx3);
%nn=max(size(yy))

idx=0;
for i=1:n
 if (td(i)>0.001)
  idx=idx+1;
  num(idx)=idx; T(idx)=tT(i); type(idx)=ttype(i); d(idx)=td(i); pcc(idx)=tpcc(i);
  sal(idx)=tsal(i); wc(idx)=twc(i);
 end
end



%nye snelag har massefylden 200kg/m3
dns=Ms./200;
if (dns>0)
 %læs parametre i temp variabel
 tnum=num; tT=T; ttype=type; tDens=Dens; td=d; tpcc=pcc; tsal=sal; twc=wc;
 %sæt første lag
 num(1)=1; T(1)=Ta; type(1)=1; Dens(1)=200; d(1)=dns; pcc(1)=0.07; sal(1)=0; wc=0;
 for i=2:n+1
   num(i)=i; T(i)=tT(i-1); type(i)=ttype(i-1); Dens(i)=tDens(i-1); d(i)=td(i-1); 
   pcc(i)=tpcc(i-1); sal(i)=tsal(i-1); wc(i)=twc(i-1);
 end
end

% n=max(size(num));
n=length(num);

end %program slut

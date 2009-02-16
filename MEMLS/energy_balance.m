function energy_balance(ifile,mfile)

eval (['load ',mfile])
eval (['pressure = ',mfile(1:6),'(:,1);']);
mm = find(pressure > 0);
eval (['tid = ',mfile(1:6),'(mm,1);']);
eval (['P = ',mfile(1:6),'(mm,2);']);
eval (['Ta = ',mfile(1:6),'(mm,3);']);
eval (['u = ',mfile(1:6),'(mm,4);']);
eval (['Sw = ',mfile(1:6),'(mm,5);']);
eval (['Lw = ',mfile(1:6),'(mm,6);']);
eval (['r = ',mfile(1:6),'(mm,7);']);
eval (['Mr = ',mfile(1:6),'(mm,8);']);

eval (['load ',ifile])
eval (['di = ',ifile(1:6),'(:,5);']);
ii = find(di > 0);
eval (['num = ',ifile(1:6),'(ii,1);']);
eval (['T = ',ifile(1:6),'(ii,2);']);
eval (['typ = ',ifile(1:6),'(ii,3);']);
eval (['Dens = ',ifile(1:6),'(ii,4);']);
eval (['d = ',ifile(1:6),'(ii,5);']);
eval (['pcc = ',ifile(1:6),'(ii,6);']);
eval (['sal = ',ifile(1:6),'(ii,7);']);
eval (['wc = ',ifile(1:6),'(ii,8);']);
eval (['rms = ',ifile(1:6),'(ii,9);']);
eval (['si = ',ifile(1:6),'(ii,10);']);

% n=max(size(num));
n=length(num);

if n == 0
   return
end %endif

tid=clock;
stid=tid(3).*86400+tid(4).*3600+tid(5).*60+tid(6);

%latent heat of fusion: vand
Lf=0.334e6;

%lufttryk hPa->Pa
P=P.*100;

snowpool=0.0;
mm=600; %tidsstep i sekunder

%?bne udfil
[fnum,err]=fopen('emipar12.temi','w');

%--------------------------
%l?kke over meteorologiske param
%--------------------------
%for mi=1:1216
%for mi=1:1210
for mi=1:1100
 %  for mi=1:9
%greenice data 239, ecwmf data 1216.

n
mi.*0.25

r(mi)=relative_humidity(Ta(mi),r(mi));

%regn eller sne
Ms=0;
if (Ta(mi) <= 273.15) Ms=Mr(mi); Mr(mi)=0; end %endif
%if (Ta(mi) <= 274.) Ms=Mr(mi); Mr(mi)=0; endif
%if (Ms<1.00) Ms=0.; end %endif
%bagatelgr?nse for sne
%if (Mr(mi)<1.00) Mr(mi)=0.; end %endif
%bagatelgr?nse for sne, resten l?gges i snowpool
if (Ms<1.00)
  snowpool=snowpool+Ms;
  Ms=0.;
else
  Ms=Ms+snowpool;
  snowpool=0.;
end %endif

Ms=Ms;

%ryd op i lag
[n,num,T,typ,Dens,d,pcc,sal,wc,rms,si]=add_layer(n,num,T,typ,Dens,d,pcc,sal,wc,rms,si,Ms,Ta(mi),u(mi));

%---------------------
%her 6 tidstep per time
%---------------------
m=1;
%while (m<3600) %antal sekunder mellem meteorologiske oplysninger
while (m<21600) %6timer

%check for afsmeltede lag
[n,num,T,typ,Dens,d,pcc,sal,wc,rms,si]=del_layer(n,num,T,typ,Dens,d,pcc,sal,wc,rms,si);
num=masucolumn(num);
T=masucolumn(T);
typ=masucolumn(typ);
Dens=masucolumn(Dens);
d=masucolumn(d);
pcc=masucolumn(pcc);
sal=masucolumn(sal);
wc=masucolumn(wc);
rms=masucolumn(rms);
si=masucolumn(si);

Cw = wc .*d .*Dens;
cs(1:n)=0;
Qc(1:n)=0;
Qs(1:n)=0;
Qtot(1:n)=0;
snowload(1:n)=0;
wcb(1:n)=0;

%kortb?lge energi i hvert lag 
Fs=shortwave_energy1(typ,Sw(mi),d,n,Dens,pcc,mm);

%varmefylde af sne/sjap/havis
cs=s_heat_salmix(sal,T,Dens);

%-----------------
%l?kke over n lag
%-----------------

for i=1:n %for n lag

%temperaturen i laget s?ttes til smeltepunktet n?r der er vand i
if (wc(i)>0.001) T(i)=meltpoint_depres(sal(i)); end %endif

%brinepocket st?rrelsen er proportional med temperaturen.
if (typ(i)==3) pcc(i)=brinepocket_pec(Vb(T(i)-273.15,sal(i))); end %endif

%fully developed depth hoar if pcc > 0.2mm then it is type 2 snow
if (pcc(i) > 0.2 && typ(i) == 1) typ(i)=2.0; end %endif

%---------------------------
%overfladebalance per tidsenhed
%---------------------------

if (i==1)

  if (typ(1)<=2) 
%    Qc(1)=snow_conductivity(Dens(1),T(1)).*dTdz(T(1),T(2),d(1)+0.5.*d(2)); 
    Qc(1)=snow_conductivity(Dens(1),T(1)).*dTdz(T(1),T(2),d(1)); 
    wcb(1)=Vb(T(1)-273.15,sal(1));
  end %endif
  if (typ(1)>2) Qc(1) = por_ice_conductivity(sal(1),T(1),Dens(1)) .* dTdz(T(1),T(2),d(1)); end %endif
  Qcc=-mm.*Qc(1);
  Qs=Fs(1);
  Qe=latent_heat(u(mi),r(mi),Ta(mi),T(1),P(mi),mm);
  Qh=sensible_heat(u(mi),Ta(mi),T(1),mm);
  Qr=precipitation_heat(Mr(mi),Ta(mi),T(1),mm);
  Ql=longwave_energy(T(1),Lw(mi),mm);

  Qtot(1)=Qcc+Qs+Qe+Qh+Qr+Ql;
  
  if (typ(1)<=2) snowload(1)=sum(0.5.*Dens(1).*d(1).*9.8); end %endif
  
  %t?r-energi
  if (T(1) < meltpoint_depres(sal(1)) & d(1)>0.001)
    cs(1)=s_heat_salmix(sal(1),T(1),Dens(1));
     T(1)=dry_temperature(T(1),Qtot(1),cs(1),Dens(1),d(1));

     d(1)=snow_depth(d(1),Dens(1),T(1),snowload(1),mm); %skal ikke n?dvendigvis opdateres ved hvert tidsstep!
   %pcc(1)=wet_snow_growth(pcc(1),wcb(1),mm);
   %if (sal(1)==0 & m>=21000) pcc(1)=depth_hoar(pcc(1),T(1),0.01.*dTdz(T(1),T(2),(d(1)+d(2))./2),Dens(1),36.*mm); endif
     Dens(1)=ddensity(Dens(1),snowload(1),mm);
    Dens(1)=rain_mass(Dens(1),Mr(mi),d(1));
  end %endif

  %v?d-energi
  if (T(1) >= meltpoint_depres(sal(1)) & d(1)>0.001)
    if (typ(1)<=2 & Cw(1)+(Qtot(1)./Lf)>=0.0) Cw(1)=Cw(1)+(Qtot(1)./Lf); T(1) = meltpoint_depres(sal(1)); end %endif
%    else T(1)=T(1)+((Cw(1)+(Qtot(1)./Lf))./(cs(1).*Dens(1).*d(1))); Cw(1)=0.0; wc(1)=0.0; endif
    if (typ(1)>2 && cs(1).*Dens(1).*d(1)~=0) T(1)=T(1)+((Cw(1)+(Qtot(1)./Lf))./(cs(1).*Dens(1).*d(1))); Cw(1)=0.0; wc(1)=0.0; end %endif
    if (typ(1)<=2) Cw(1)=Cw(1)+Mr(mi); end %endif
    if (typ(2)<=2) Cw(1)=Cw(1)-water_flux(wc(1),mm,pec2d(pcc(1)),Dens(1)); end %endif
    if (d(1)>0.001 && Dens(1)>0) wc(1)=Cw(1)./(Dens(1).*d(1)); end %endif
    pcc(1)=wet_snow_growth(pcc(1),wc(1),mm);
    d(1)=snow_depth(d(1),Dens(1),T(1),snowload(1),mm);
    Dens(1)=ddensity(Dens(1),snowload(1),mm);
    [Dens(1),d(1)]=snow_melt(Dens(1),d(1),Qtot(1),Mr(mi),-water_flux(wc(1),mm,pec2d(pcc(1)),Dens(1)));
    end %endif

end %endif % if surface layer 

%----------------
%mellemliggende lag
%----------------

if (i>1 & i<n)
  if (typ(i)<=2) 
%    Qc(i)=snow_conductivity(Dens(i),T(i)).*dTdz(T(i-1),T(i+1),d(i)+0.5.*d(i-1)+0.5.*d(i+1)); 
%    Qc(i)=snow_conductivity(Dens(i),T(i)).*dTdz(T(i-1),T(i+1),d(i)); 
    Qctop=snow_conductivity(Dens(i),T(i)).*dTdz(T(i-1),T(i),d(i));
    Qcbot=snow_conductivity(Dens(i),T(i)).*dTdz(T(i+1),T(i),d(i));
    wcb(i)=Vb(T(i)-273.15,sal(i));
  end %endif
  if (typ(i)>2) 
%    Qc(i)=por_ice_conductivity(sal(i),T(i),Dens(i)).*dTdz(T(i-1),T(i+1),d(i)+0.5.*d(i-1)+0.5.*d(i+1)); 
%    Qc(i)=por_ice_conductivity(sal(i),T(i),Dens(i)).*dTdz(T(i-1),T(i+1),d(i)); 
    Qctop=por_ice_conductivity(sal(i),T(i),Dens(i)).*dTdz(T(i-1),T(i),d(i));
    Qcbot=por_ice_conductivity(sal(i),T(i),Dens(i)).*dTdz(T(i+1),T(i),d(i));
  end %endif
  Qs=Fs(i);
  Qcc=mm.*(Qctop+Qcbot);
%  Qcc=mm.*Qc(i);
  Qtot(i)=Qs+Qcc;

  snowload(i)=sum(Dens(1:i-1).*d(1:i-1).*9.8);

  %t?r-energi
  if (T(i) < meltpoint_depres(sal(i)) & d(i)>0.001)
    cs(i)=s_heat_salmix(sal(i),T(i),Dens(i));
     T(i)=dry_temperature(T(i),Qtot(i),cs(i),Dens(i),d(i));
    if (typ(i)<=2) 
      d(i)=snow_depth(d(i),Dens(i),T(i),snowload(i),mm);
      Dens(i)=ddensity(Dens(i),snowload(i),mm);
     % pcc(i)=wet_snow_growth(pcc(i),wcb(i),mm); 
      if (sal(i)==0 & m>=21000 & i>=3 & d(i)>0.005)
        pcc(i)=depth_hoar(pcc(i),T(i),0.01.*dTdz(T(i-1),T(i+1),d(i)+0.5.*d(i-1)+0.5.*d(i+1)),Dens(i),36.*mm);
      end %endif
    end %endif
  end %endif

  %v?d-energi
  if (T(i) >= meltpoint_depres(sal(i)) & d(i)>0.001)
    if (typ(i)<=2 & Cw(i)+(Qtot(i)./Lf)>=0.0) Cw(i)=Cw(i)+(Qtot(i)./Lf); T(i)=meltpoint_depres(sal(i));
    else T(i)=T(i)+((Cw(i)+(Qtot(i)./Lf))./(cs(i).*Dens(i).*d(i))); Cw(i)=0; wc(i)=0.0; end
    if (typ(i+1)<=2) Cw(i)=Cw(i)-water_flux(wc(i),mm,pec2d(pcc(i)),Dens(i))+water_flux(wc(i-1),mm,pec2d(pcc(i-1)),Dens(i-1)); end %endif
    if (d(i)>0.001) wc(i)=Cw(i)./(Dens(i).*d(i)); end %endif
    if (typ(i)<=2) 
     pcc(i)=wet_snow_growth(pcc(i),wc(i),mm);
     d(i)=snow_depth(d(i),Dens(i),T(i),snowload(i),mm);
     Dens(i)=ddensity(Dens(i),snowload(i),mm);
     [Dens(i),d(i)]=snow_melt(Dens(i),d(i),Qtot(i),0,water_flux(wc(i-1),mm,pec2d(pcc(i-1)),Dens(i-1))-water_flux(wc(i),mm,pec2d(pcc(i-1)),Dens(i-1)));
    end %endif
  end %endif
end %endif

%----------------
%bundlag
%----------------

if (i==n)
  Qctop=por_ice_conductivity(sal(n),T(n),Dens(n)).*dTdz(T(n-1),T(n),(d(n-1)+d(n))./2);
  Qcthrough=por_ice_conductivity(sal(n),T(n),Dens(n)).*dTdz(T(n-1),meltpoint_depres(32),d(n));
  Qw=1.988; %maykut&untersteiner,1971 [W/m2] assumed constant throughout year
  Qc(n)=por_ice_conductivity(sal(n),T(n),Dens(n)).*dTdz(T(n-1),meltpoint_depres(32),d(n));
  Qcc=mm.*(Qcthrough+Qw);
  Qs=Fs(n);
  Qtot(n)=Qcc+Qs;
  cs(n)=s_heat_salmix(sal(n),T(n),Dens(n));

  dfor=d(n);
  d(n)=d(n)-(Qtot(n)./(Lf.*Dens(n)));
  if (Qtot(n)<0)
    [salbot,dbot]=newlayer_sal(Qtot(n),Lf,Dens(n));
    sal(n)=(dfor.*sal(n)+dbot.*salbot)./(dfor+dbot);
  end
%  if (d(n)>0.001) T(n)=T(n)+((Qtot(n))./(cs(n).*Dens(n).*d(n))); endif
     T(n)=dry_temperature(T(n),Qtot(n),cs(n),Dens(n),d(n));
end %endif

end %endfor %for n lag

m=m+mm;

end %endwhile %while m iterationer over T

di=flipud(d);
Ti=flipud(T);
Wi=flipud(wc);
roi=flipud(Dens);
pci=flipud(pcc);
sali=flipud(sal);
typei=flipud(typ);
sii=flipud(si);

[md]=memls_mod(num,di,Ti,Wi,roi,pci,sali,typei,sii);

% check processing time
tid=clock;
ptid=tid(3).*86400+tid(4).*3600+tid(5).*60+tid(6)-stid

[sfb,dlay]=freeboard(d,Dens);
[snowt,icet]=thickness(n,d,typ);
[adens,apcc,asT,asal,aiT,isT,pcc1]=snowavg(n,typ,d,Dens,pcc,T,sal);
%typ
T

page_output_immediately=1;
%fflush(fnum);

data=real([mi;Ta(mi);Ms+Mr(mi);sfb;snowt;icet;d(1);Dens(1);T(1);pcc1;adens;apcc;asT;asal;aiT;isT;md(1);md(2);md(3);md(4);md(5);md(6);md(7);md(8);md(9);md(10);md(11);md(12);md(13);md(14);md(15);md(16);md(17);md(18);md(19);md(20)]);
fprintf(fnum,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',data');
%data'

[fnum3,err]=fopen('dump12.temi','a');
page_output_immediately=1;
%fflush(fnum3);
data3=[num,T,typ,Dens,d,pcc,sal,wc];
fprintf(fnum3,'%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',0.0,mi,Ta(mi),adens,apcc,snowt,sfb,T(1));
if (rem(mi,10)==0)
fprintf(fnum3,'%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',data3');
end %endif
fclose(fnum3)

end %endfor %for met obs

fclose(fnum)

end %program slut

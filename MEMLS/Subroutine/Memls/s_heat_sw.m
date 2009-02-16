function shsw=s_heat_sw(s,t,p)

%salinity,s[psu],temperature,t[C],pressure,p[bar](atmospheric pressure=0)

%computes the specific heat of sea water
%Fofonoff, practice of oceanography,JGR 90,c2,3332-3342,1985
%Fofonoff & Millard, UNESCO, 1983.

sr=sqrt(s);

%cp0
a=(-1.38385e-3.*t+0.1072763).*t-7.643575;
b=(5.148e-5.*t-4.07718e-3).*t+0.1770383;
c=(((2.093236e-5.*t-2.654387e-3).*t+0.1412855).*t-3.720283).*t+4217.4;
cp0=(b.*sr+a).*s+c;

%cp1
a=(((1.7168e-8.*t+2.0357e-6).*t-3.13885e-4).*t+1.45747e-2).*t-0.49592;
b=(((2.2956e-11.*t-4.0027e-9).*t+2.87533e-7).*t-1.08645e-5).*t+2.4931E-4;
c=((6.136e-13.*t-6.5637e-11).*t+2.6380e-9).*t-5.422e-8;
cp1=((c.*p+b).*p+a).*p;

a=(((-2.179e-10.*t+2.5941e-8).*t+9.802e-7).*t-1.28315e-4).*t+4.9247e-3;
b=(3.122e-8.*t-1.517e-6).*t-1.2331e-4;
a=(a+b.*sr).*s;
b=((1.8448e-11.*t-2.3905e-9).*t+1.17054e-7).*t-2.9558e-6;
b=(b+9.971e-8.*sr).*s;
c=(3.513e-13.*t-1.7682e-11).*t+5.540e-10;
c=(c-1.43e-12.*t.*sr).*s;
cp2=((c.*p+b).*p+a).*p;

shsw=cp0+cp1+cp2;

end %endfunction

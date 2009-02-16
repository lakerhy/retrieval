function extinction=shortwave_extinction(type,n,rho,pec)

%something maybe found in boren & barkstrom, 1972,grenfell&maykut, 1977,maykut, 1986; grenfell& perowich, 1984

%formulation of extinction coefficien Jordan et al. 1999 JGR eq.12
%beta=0.003795.*rho./sqrt(d);

%grain diameter
gs=d2do(pec2d(pec),rho);

extinction(1:n)=0;
lookup=[30,20,2,2.5];

for i=1:n
 if (type(i)<=2)
  extinction(i)=0.003795.*rho(i)./sqrt(gs(i));
 else
  extinction(i)=lookup(type(i));
 end %endif
end %endfor
end

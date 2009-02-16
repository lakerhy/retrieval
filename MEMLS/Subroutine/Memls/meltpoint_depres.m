function m=meltpoint_depres(S)
%Tf=Kfm Granberg 1998 eq. 5, smeltepunktssænkning af sne i kontakt med brine...?
%molvægten af NaCl: 58.45g/mol
%vægt af cm3 is: 0.92g
%frysepunktssænkningskonstant: 1.86K/mole
m=273.15-0.05411.*S; %Maykut (1986) p. 433
%m=273.15-(1.86.*0.92.*S./58.45);
%m=273.15-0.0293.*S; %Granberg angiver dette for sne

end

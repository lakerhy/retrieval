function albedo=shortwave_albedo(type)

%something maybe found in boren & barkstrom, 1972,maykut, 1986; grenfell& perowich, 1984

albedo=0;

lookup=[0.9,0.8,0.75,0.75];
albedo=lookup(type);

end

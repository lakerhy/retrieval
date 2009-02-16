function [albedo,extinction]=shortwave_spectral1(type,n)

%something maybe found in boren & barkstrom, 1972,maykut, 1986; grenfell& perowich, 1984

albedo(1:n)=0;
extinction(1:n)=0;

for i=1:n
%albedo
if type(i)==1 albedo(i)=0.9; end
if type(i)==2 albedo(i)=0.8; end
if type(i)==3 albedo(i)=0.75; end
if type(i)==4 albedo(i)=0.75; end

%extinction
if type(i)==1 extinction(i)=30; end
if type(i)==2 extinction(i)=20; end
if type(i)==3 extinction(i)=2; end
if type(i)==4 extinction(i)=2.5; end

end
%call this func: [a,b,c]=shortwave_spectral1(type,n)
end

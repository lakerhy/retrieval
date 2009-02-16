function [sal,d]=newlayer_sal(Qtot,Lf,Dens)
sal=0.0;
d=-(Qtot./(Lf.*Dens));
grate=d/6.0; %v?kstrate i cm/s, bestemmer saliniteten
if (d>0.0)
%Nakawo & Sinha, 1981, figur 13.
sal=32.0.*0.12./(0.12+0.88.*exp(-4.2e4.*grate));
end %endif

end %endfunction

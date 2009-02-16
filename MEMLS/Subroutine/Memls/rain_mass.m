function dens=rain_mass(dens,mr,dz)
 
if (dz<0.001) 
    return 
end

dens=dens+mr./dz;

if (dens>926) 
    dens=926; 
end

end

function lw=longwave_energy(Ts,lw_i,time)
%computes the net longwave energy at the surface brun et al. 1989
sbc=5.67051e-8; %stephan-boltzmans constant
epsilon_s=0.98; %surface emissivity

lw=time.*((epsilon_s.*lw_i)-((epsilon_s.*sbc).*(Ts.^4)));

end

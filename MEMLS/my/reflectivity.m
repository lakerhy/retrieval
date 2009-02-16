function[] = reflectivity(eps1,eps2,theta1,theta2,freq,dei)
  gamma = abscoeff(real(eps1),imag(eps1),0,freq,0)
  ri =0
  ti = exp(gamma*dei.*(-1));
  sih(1) = 
  roa_h = ri + sih(1)*ti^2/(1-ri*sih(1));
   
      rob_h = ri + sih(2)*ti^2/(1-ri*sih(2));
      
      rsa_h(n,i) = sih(2)+(1-sih(2))^2*roa_h./(1-sih(2)*ri);
      
      rsb_h(n,i) = sih(1)+(1-sih(1))^2*roa_h./(1-sih(1)*ri);
  

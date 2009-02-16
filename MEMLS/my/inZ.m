function[Zin] = inZ (Z,R,gamma,tei)
  %Z transmission line intrinsic impedance
    % R reflection coefficient
    % tei
    
    Zin = Z.*((1+R.*exp(-2*i*gamma*tei))./(1-R.*exp(-2*i*gamma*tei)));

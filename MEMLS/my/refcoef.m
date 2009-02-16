function[R] = refcoef (Z1,Z2)
  % calculate the reflection coefficients
    %Z1 is the intrinsic impedance of transmission line and Z2 is the
    %intrinsic coefficient of load
    
    R = zeros(1,2);
    R(1)= (Z1(1)-Z2(1))./(Z1(1)+Z2(1));%V
    R(2) = (Z2(2)-Z1(2))./(Z2(2)+Z1(2)); %H

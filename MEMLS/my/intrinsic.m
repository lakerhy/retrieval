function[Z] = intrinsic (epsi,theta)
  
  % Calculate the intrinsic impedance
    % epsi: dielectric  constant of the media
    % theta: incidental angle in the media [rad]
    Z = zeros(1,2);
Z(1) = 1/sqrt(epsi)*cos(theta);% V;
Z(2) = 1/sqrt(epsi)*sec(theta);% H;

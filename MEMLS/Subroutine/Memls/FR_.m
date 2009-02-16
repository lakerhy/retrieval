%computes the Fresnell reflection coefficient
function r=FR_(ia,ta,e1,e2)
nn1=1. ./sqrt(real(e1));
nn2=1. ./sqrt(real(e2));
r=(nn2 .*cos(ia)-nn1 .*cos(ta)) ./(nn2 .*cos(ia)+nn1 .*cos(ta));

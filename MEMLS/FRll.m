%computes the Fresnell reflection coefficient
function r=FRll(ia,ta,e1,e2)
nn1=1. ./sqrt(real(e1));
nn2=1. ./sqrt(real(e2));
r=(nn1 .*cos(ia)-nn2 .*cos(ta)) ./(nn1 .*cos(ia)+nn2 .*cos(ta));

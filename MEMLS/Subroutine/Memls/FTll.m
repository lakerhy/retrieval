function Tll=FTll(ia,ta,e1,e2)
%computes the transmission coefficient (VV) using the angle and dielectrics
 nn1=1.0 ./sqrt(real(e1));
 nn2=1.0 ./sqrt(real(e2));
 Tll=(2.0 .*nn2 .*cos(ia)) ./(nn2 .*cos(ta)+nn1 .*cos(ia));

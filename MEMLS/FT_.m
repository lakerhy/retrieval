function T_=FT_(ia,ta,e1,e2)
%computes the Transmission coefficient (HH), using angles and dielectrics
 nn1=1.0 ./sqrt(real(e1));
 nn2=1.0 ./sqrt(real(e2));
 T_=(2.0 .*nn2 .*cos(ia)) ./(nn2 .*cos(ta)+nn1 .*cos(ta));

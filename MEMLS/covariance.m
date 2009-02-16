function[cov] = covariance(dev_corr)
%  compute the covariance matrix from deviation and correlation
  
  
  
  dev = dev_corr(1,:);
  correlation = dev_corr(2:8,:);
 
  for i = 1:7 
    for j = 1:i
    cov(i,j) = correlation(i,j)*dev(i)*dev(j);
    cov(j,i) =  cov(i,j);
    
    end 
  end
  
  
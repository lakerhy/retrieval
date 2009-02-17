function[] = dayretrieval (fname)
  
  cd ../MEMLS/Data
  temp = load(fname);
  cd ../../inversion
  data = temp(:,3:12);
  [a,b] = size(data)
  for i=1:a
p_est(i,:)= inversion(data(i,:));
  end
  
  plot(1:155,p_est(:,2));
  

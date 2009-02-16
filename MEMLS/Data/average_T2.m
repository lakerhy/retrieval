function [output,count] = average_T2(fname)
  
  

  A = load(fname);
  count = 0;
  index1 = find(A(:,2)<3);
  output = zeros(1,10);
  % brightness temperature of AMSR oberservation rows:Jan to Sep cols: frequency bands 6V, 6H, 10V, 10H, 19V, 19H, 23V, 23H, 37V, 37H
  %output = zeros(9,10);
  [rows, cols] = size(A); 
  for i=1:length(index1)
    
    output = output + A(i,6:15);
    
    count = count+1;
  end
  [rows,cols]= size(A);
  for i = 1:rows
    if((A(i,2)==3)&&(A(i,3)<16))
      output = output + A(i,6:15);
      count = count +1;
    end
  end
  
  output = output/count;
  
  
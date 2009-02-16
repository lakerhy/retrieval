function [output,count] = average_T(fname)
  
  

  A = load(fname);
  count = 0;
  B = A(:,1)/1e4;
  B = floor(B)-20040116;
  index = find(B(:,1)<0);
  output = zeros(10)
  % brightness temperature of AMSR oberservation rows:Jan to Sep cols: frequency bands 6V, 6H, 10V, 10H, 19V, 19H, 23V, 23H, 37V, 37H
  %output = zeros(9,10);
  [rows, cols] = size(A); 
  for i=1:length(index)
              output = output + A(i,2:cols-3);
        count = count+1;
   end
  
     output = output/count;
 
  
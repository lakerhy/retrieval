function [output,count] = average_T(fname,type)

  A = load(fname);
  count = 0;
  output = zeros(1,10);
  index = 0;
  start = 0;

  if type == 1; %% MY 

%     B = A(:,1)/1e4;
%     B = floor(B)-20040116;
%     index = find(B(:,1)<0);
%     % brightness temperature of AMSR oberservation rows:Jan to Sep cols: frequency bands 6V, 6H, 10V, 10H, 19V, 19H, 23V, 23H, 37V, 37H
%     %output = zeros(9,10);
    
    
%     start = 2;
%       for i=1:length(index)
%     output = output + A(i,start:start+9);
%     count = count+1;
%       end
%       output = output/count;
      
     
    [row,col] = size(A);
    
    start=2;
 
    A_sum = sum(A);
    
    output = A_sum(:,start:start+9)/row;
    
  end

  

  if type ==2; %% FY
    
    %B1 = find(A(:,2)==1);
    %B2 = find(A(1:length(B1),3)<=15);
    
    [row,col] = size(A);
    
    start=6;
 
    A_sum = sum(A);
    
    output = A_sum(:,start:start+9)/row;
    end
  
 
  

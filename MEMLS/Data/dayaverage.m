% Write the day average Tb into file dailyTb
function[count] = dayaverage (filename)

  data = load(filename);
  [a,b]=size(data);

  count = zeros(5,31);
  Tb = zeros(10,31*5);
  for month = 1:5;
    for day = 1:31;
      for i=1:a
        if (month==data(i,2)&& day==data(i,3))  
          count(month,day) = count(month,day)+1;
        end
      end
    end
  end
  mark = 1;
  for month = 1:5
    for day = 1:31
      if count(month,day)>0
      buffer = data(mark:mark+count(month,day),6:15);
      %Tb(month,day,:) = mean(buffer);
      Tb(:,31*(month-1)+day)=mean(buffer);
      mark = mark+count(month,day);
      end
    end
  end
  fmonth = zeros(1,month*day);
    fday = zeros(1,month*day);
  fmonth(1,1:31)        = 1*ones(1:31,1);
  fmonth(1,31*1+1:31*2) = 2*ones(1:31,1);
  fmonth(1,31*2+1:31*3) = 3*ones(1:31,1);
  fmonth(1,31*3+1:31*4) = 4*ones(1:31,1);
  fmonth(1,31*4+1:31*5) = 5*ones(1:31,1);
  fday(1,1:31) = [1:31]
  fday(1,31*1+1:31*2) = [1:31];
  fday(1,31*2+1:31*3) = [1:31];
  fday(1,31*3+1:31*4) = [1:31];
  fday(1,31*4+1:31*5) = [1:31];
  foutbuffer= [fmonth;fday;Tb]
  fid = fopen('dailyTb', 'w');
 
fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',foutbuffer);
fclose(fid)
  
    

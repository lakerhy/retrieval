function[] = dayretrieval2
  fname1='area.58.2005.1-6';
fname2 = 'area.33.2005.1-6';
  cd ../MEMLS/Data
  tempFY = load(fname1);
  tempMY = load(fname2);
  cd ../../inversion2
  dataFY = tempFY(:,6:15);
  dataMY = tempMY(18207:28409,6:15);
  [aFY,bFY] = size(dataFY);
  [aMY,bMY] = size(dataMY);

  day1=aFY;
  day2=aMY;
  for i=1:day1
   [p_est_FY(i,:),std_S_FY(i,:),std_Sp_FY(i,:),DeltaTb_FY(i,:)]=  
inversion(dat$

  end


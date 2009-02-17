function[] = dayretrieval 
  
  cd ../MEMLS/Data
  tempFY = load('58.2005.1.1.7');
  tempMY = load('33.2005.1.2.7');
  cd ../../inversion2
%   dataFY = tempFY(:,3:12);
%   dataMY = tempMY(:,3:12);
 dataFY = tempFY(:,6:15);
 dataMY = tempMY(:,6:15);
 [aFY,bFY] = size(dataFY);
  [aMY,bMY] = size(dataMY);
  
  day1=aFY;
  day2=aMY;
  for i=1:day1
   [p_est_FY(i,:),std_S_FY(i,:),std_Sp_FY(i,:),DeltaTb_FY(i,:)]=  inversion(dataFY(i,:))
 
  end
  
    fid =fopen('P_est_FY','w');
  fprintf(fid,'%d,%d,%d,%d,%d\n',p_est_FY');
  fclose(fid);
  
  fid=fopen('std_S_FY','w');
  fprintf(fid,'%d,%d,%d,%d,%d\n',std_S_FY');
  fclose(fid);

  
  fid=fopen('DeltaTb_FY','w');
  fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',DeltaTb_FY');
  fclose(fid);


  
  for i=1:day2
    [p_est_MY(i,:),std_S_MY(i,:),std_Sp_MY(i,:),DeltaTb_MY(i,:)]= inversion(dataMY(i,:))
  end
  
  fid =fopen('P_est_MY','w');
  fprintf(fid,'%d,%d,%d,%d,%d\n',p_est_MY');
  fclose(fid);
  
  
  fid=fopen('std_S_MY','w');
  fprintf(fid,'%d,%d,%d,%d,%d\n',std_S_MY');
  fclose(fid);
  
  fid=fopen('DeltaTb_MY','w');
  fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',DeltaTb_MY');
  fclose(fid);


  

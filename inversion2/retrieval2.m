function[] = retrieval2(type)

  if type == 1 
    % Inversion evaluation
    %p=[250,18,0.1,320,0]'; 
    p=[255,20,0.14,320,0.45,0.26]';
    Tb = fw(p);
    
  end

  if type == 2
    % FY data 58s    
 %   [Tb_AMSR]=amsr('area.46.2005.sort');
      [Tb_AMSR]=amsr('area.58.2005.sort');
    Tb = Tb_AMSR;
  end
  
  if type ==3
    % MY data 33s    
%    [Tb_AMSR]=amsr('area.33.2005.sort');
       [Tb_AMSR]=amsr('area.33.2005.sort');
    Tb = Tb_AMSR;
  end

  [p_est,S_std,Sp_std]=inversion(Tb)


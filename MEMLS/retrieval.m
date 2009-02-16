
function[p] = retrieval(type)

  inv = 1
  if inv==1 % using inv_funtest (via fw_fun2)
    
    if type==1%% MY
      cd Data;
      output = average_T('area.33s.sort',1);
      cd ..;
      
      [p,S]=inv_funtest(output,type);
      %[p,S]=inv_fun(output);
    end
    
    if type==2 %%FY
      cd Data;
      output = average_T('area.58s.sort',2);
      cd ..;

      [p,S]=inv_funtest(output,type);
    end
    

    S
    
  else % using inv_fun (via fw_fun3)
    
    if type==1%% MY
      cd Data;
      output = average_T('area.33s.sort',1);
      cd ..;
      
      [p,S]=inv_funtest(output,type);
      %[p,S]=inv_fun(output);
    end
    
    if type==2 %%FY
      cd Data;
      output = average_T('area.58s.sort',2);
      cd ..;

      [p,S]=inv_fun(output,type);
    end
    
s    Tb_plot_fun(type,p,output);
    
  end
  
  
  

% ******** Usage: Compute  Emissivity of MEMLS model
% ******** Input ********************************** 
% memls_profile: memls snow-ice profile file name
% type: 3/4 FY/MY
% ******** Onput **********************************
% emissivity 5x2

function [emissivity] = emi(par,type)
  
  Tb_100K = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
                             par,type,100);
  
  Tb_0K = icemain('Freq-Memls-in.txt','Angl-Memls-in.txt', ...
                             par,type);
  
  emissivity = 1 - (Tb_100K-Tb_0K)/100;
  
%   fid = fopen(output_file,'wt');
%  fprintf(fid,'%4.15f \n', output);
%  fclose(fid);

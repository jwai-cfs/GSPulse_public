function s = load_json_dict(fn)
% Description: 
%  load a .json file into a matlab struct
% 
% Inputs: 
%  fn - the .json file to load
%
% Outputs:
%  s - struct containing data from .json file

fid = fopen(fn); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
s = jsondecode(str);
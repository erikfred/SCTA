function struc = merge_oneElementStructure(struc1, struc2, avoid)
% Function to merge common VECTOR fields of a 1 element structure
%
% Usage
%   struc = merge_oneElementStructure(struc1, struc2)
% 
% Inputs
%   STRUC1 - One element structure
%   STRUC2 - Another one element structure
%   AVOID  - Names of field in STRUC1 to keep and not to merge
%
% Outputs
%   STRUC - Merged structure

if nargin<3
  avoid = '';
end

if nargin~=3 && nargin~=2
  error('merge_oneElementStructure requires 2-3 arguments')
end

f1 = fieldnames(struc1);
f2 = fieldnames(struc2);
struc = struc1;

for i = 1:length(f1)
  if ~any(strcmp(f1{i},avoid))
    
    if any(strcmp(f1{i},f2))

      eval(['dim1 = size(struc1.' f1{i} ');'])
      eval(['dim2 = size(struc1.' f1{i} ');'])

      if dim1(1)==1 && dim2(1)==1;
        eval(['struc.' f1{i} ' = [ struc1.' f1{i} ' struc2.' f1{i} '];'])
      else dim2(1)==1 && dim2(1)==1;
        eval(['struc.' f1{i} ' = [ struc1.' f1{i} '; struc2.' f1{i} '];'])
      end
    end
  end
end

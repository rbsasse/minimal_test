function [input] = util_expand(params)
%UTIL_EXPAND Generate all combinations of parameters
%
% Given an input structure with names and numeric values
% this function returns a structure with all combinations
% which may be used for running simulations over multiple parameters
%
% PARAMETERS:
% - params: structure of named numeric parameters
%
% RETURNS: structure with named parameters
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

nm = fieldnames(params);
params = struct2cell(params);
input = cell(1,numel(params));
[input{:}] = ndgrid(params{:});
n = length(input);
for k=1:n
    input{k} = input{k}(:);
end
input = cell2struct(input, nm, 2);
end

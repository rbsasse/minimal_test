function [out, hash] = util_filename(input, prefix, postfix)
%UTIL_FILENAME Automatic filename
%
% Generate a filename based on a data structure with numeric values
%
% PARAMETERS:
% - input structure of named parameters
% - prefix [optional]
% - postfix [optional]
%
% RETURNS: string, and optional md5 hash for a more concise filename
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

if(nargin < 3)
    postfix='';
end
if(nargin < 2)
    prefix='';
end

fields = fieldnames(input);
n = numel(fields);
fname = '';
for k = 1:n
    par = fields{k};
    val = input.(par);
    if(ischar(val))
        toadd = sprintf('%s=%s', par, val);
    elseif (islogical(val) && val) % boolean true
        toadd = sprintf('%s', par);
    elseif (islogical(val) && ~val) % boolean false
        toadd = sprintf('not-%s', par);
    elseif (isnumeric(val) && length(val)==1)
        toadd = sprintf('%s=%g', par, val);
    elseif (isnumeric(val) && length(val)>=1)
        toadd = sprintf('%s=%g-%g', par, min(val), max(val));
    end
    fname = [fname, toadd, '_'];
end
out = [prefix, fname, postfix, 'auto.mat'];

% if we request a md5 hash, get it from a system call (not on windows)
if(nargout > 1)
    out = [prefix, fname, postfix];
    hash = util_md5_external(out);
end

end

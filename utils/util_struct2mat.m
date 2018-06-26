function [ m ] = util_struct2mat( s, direction )
%UTIL_STRUCT2MAT Convert structure to matrix
%
% Given a structure of matrices with compatible sizes, binds them in a
% matrix format.
%
% PARAMETERS:
% - s structure of matrices
% - direction [optional] bind vertically or horizontally
%
% RETURNS: matrix
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

if(nargin < 2)
    direction = 'v';
end

tmp = struct2cell(s);
switch direction
    case 'h'
        m = horzcat(tmp{:});
    case 'v'
        m = vertcat(tmp{:});
end
end

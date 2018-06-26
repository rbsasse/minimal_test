function [ x0 ] = inversion_cg_multidir( A, y, tol, nmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    tol=1e-4;
    nmax=10;
end

B = ctranspose(A); % Hermitian transpose.
% Note: diag blocks not symmetric, so no easier way
z = B * y;

% initial steps
x0 = y; % allocate result with first born approx (equiv. preconditioner)

p_i = z - B*A*x0; gi = p_i;
x0 = x0;
wi = A * x0;
vi = A * p_i;

n = 0;
rel_error = Inf; % large value
while ((n < nmax) && (rel_error > tol))
    
    alphai = sum(gi .* conj(gi)) ./ sum(vi .* conj(vi));
    xip1 = x0 + bsxfun(@times, p_i, alphai);
    wip1 = wi + bsxfun(@times, vi, alphai);
    gip1 = z - B * wip1;
    betai = sum(gip1 .* conj(gip1)) ./ sum(gi .* conj(gi));
    
    pip1 = gip1 + bsxfun(@times, p_i, betai);
    vip1 = A * pip1; 
    rel_error = max(norm(xip1 - x0) ./ norm(xip1)) ;
    
    if(rel_error < tol)
        x0 = xip1;
        break; % skip cost of next line
    end
    % updating
    x0 = xip1;
    p_i = pip1;
    wi = A * x0;
    vi = vip1;
    gi = gip1;
    n = n+1;
end % while

end


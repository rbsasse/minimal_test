function [ x0 ] = inversion_cg( A, y, tol, nmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    tol=1e-4;
    nmax=10;
end

[N, nr] = size(y); % number of dipoles (x3), and angles to solve for
B = ctranspose(A); % Hermitian transpose.
% Note: diag blocks not symmetric, so no easier way

% initial steps
x0 = y; % allocate result with first born approx (equiv. preconditioner)
% z = B * y;
% p0 = z - B*A*x0;
% w0 = A * x0;
% v0 = A * p0;

% initialise some complex vectors
% y_i=zeros(N,1)+0i; 

% loop over RHS cols
% (might be avoided?)
for (ii=1:nr)
    
    % initial step
    x_i = x0(:,ii);
    y_i = y(:,ii);
    z_i = B * y_i;
    g_i = z_i - B * A * x_i;
    p_i = g_i;
    w_i = A * x_i;
    v_i = A * p_i;
    
    n = 0;
    rel_error = 1e4; % large value
    while ((n < nmax) && (rel_error > tol))
        
        alpha_i = dot(g_i,g_i) / dot(v_i, v_i);
        x_ip1 = x_i + alpha_i * p_i;
        w_ip1 = w_i + alpha_i * v_i;
        g_ip1 = z_i - B * w_ip1;
        beta_i = dot(g_ip1, g_ip1) / dot(g_i, g_i);
        
        p_ip1 = g_ip1 + beta_i * p_i;
        % v_ip1 = A * g_ip1 + beta_i * v_i;
        v_ip1 = A * p_ip1; % more accurate according to Draine
        rel_error = norm(x_ip1 - x_i) / norm(x_ip1) ;
        
        if(rel_error < tol) 
            x_i = x_ip1;
            break; % skip cost of next line
        end 
        % updating
        x_i = x_ip1;
        p_i = p_ip1;
        w_i = A * x_i;
        v_i = v_ip1;
        g_i = g_ip1;
        n = n+1;
    end % while
    x0(:,ii) = x_i;
end % loop

end


function  E = solve_iterative(A,Einc)

E = Einc;

Na = size(Einc,2);
tol=1e-3;maxit=20;

for ii=1:Na
    [x,~] = bicgstab(A,Einc(:,ii),tol,maxit,[],[],Einc(:,ii));
    E(:,ii) = x;
end
end
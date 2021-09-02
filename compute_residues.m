
function res = compute_residues(A,B,C)

% This program computes the matrix residues of the partial
% fraction expansion of H(s) = C(sI- A)^-1 B.
%
% A = matrix of dimensions n x n
% B = matrix of dimensions n x p
% C = matrix of dimensions m x n
% V = Vandermonde matrix
% root = a vector containing the nonrepeated eigenvalues of A
% porder = a vector representing the order of each corresponding eigenvalue.
% res = the matrix residues of the expansion.
%

n = length(A);
dimB = size(B);
p = dimB(2);
dimC = size(C);
m = dimC(1);

e = eig(A);
eps = 1e-6; sigma = 0; i = 1;
while i <= n
    sigma = sigma + 1;
    root(sigma) = e(i);
    porder(sigma) = 1;
    for j = i + 1 : n
        if abs(real(e(i)) - real(e(j))) < eps  && abs(imag(e(i)) - imag(e(j))) < eps
            dummy = e(j); e(j) = e(i); e(i) = dummy;
            porder(sigma) = porder(sigma) +1;
            i=i+1;
        end
    end
    i=i+1;
end

%
% Construct the Vandermonde matrix from the eigenvalues of A
%
r = 0;
for k = 1: sigma
    r = r+1; van(1,r) = 1.0;
    for i = 2 : n
        van(i, r) = van(i-1, r)*root(k);
    end
    j = 2;
    while porder(k) >= j
        r=r+1;
        for i = 1 : n
            if i< j
                van(i, r) = 0;
            elseif i == j
                van(i,r) = 1;
            else
                van(i,r) = van(i-1,r-1)*(i-1)/(j-1);
            end
        end
        j=j+l;
    end
end

% Calculate the residues using the new formula
%
b(1:n,1:p) = B;
for i=1:n-1
    b(i*n+1:(i+1)*n,:) = A*b((i-1)*n+1:i*n,:);
end
res = kron(inv(van),C)*b;
disp("the residues of the expansion are =")
r = 0;
for k = 1:sigma
    fprintf("\nthe residues for the eigenvalue %g %gi \n\n",real(root(k)),imag(root(k)))
    for i = 1:porder(k)
        r=r+1;
        fprintf("factor = %2.0f\n\n",i)
        disp(res((r-1)*m+1:r*m,:))
    end
end
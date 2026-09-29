function [x, w, D] = sem_gll(p)
%SEM_GLL Legendre-Gauss-Lobatto points, weights and nodal derivative matrix.
% Adapted from MeshClass.lglnodes in this repository.
validateattributes(p, {'numeric'}, {'scalar','integer','>=',1});
if p == 1
    x = [-1; 1];
    w = [1; 1];
else
    x = -cos(pi*(0:p)'/p);
    for iteration = 1:100
        P0 = ones(size(x));
        P1 = x;
        for k = 2:p
            P2 = ((2*k-1)*x.*P1-(k-1)*P0)/k;
            P0 = P1;
            P1 = P2;
        end
        dx = (x.*P1-P0)./((p+1)*P1);
        x = x-dx;
        if max(abs(dx)) < 5e-15
            break
        end
    end
    w = 2./(p*(p+1)*P1.^2);
end
% Barycentric differentiation; D(i,j) = l_j'(x_i).
b = ones(p+1,1);
for j = 1:p+1
    b(j) = 1/prod(x(j)-x([1:j-1,j+1:end]));
end
D = zeros(p+1);
for i = 1:p+1
    for j = 1:p+1
        if i ~= j
            D(i,j) = b(j)/(b(i)*(x(i)-x(j)));
        end
    end
    D(i,i) = -sum(D(i,:));
end
end

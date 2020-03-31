function A_xi = S2C_matrix(n,xi)
% This function A_xi = S2C_matrix(n,xi) generate matrix A_xi which is in a
% star-to-complete demand pattern with n locations and xi. A_xi =
% (1-xi)*A_S + xi*A_C.
A_S = zeros(n);
A_S(2:n,1) = 1;
A_S(1,2:n) = 1/(n-1);
A_C = 1/(n-1)*(ones(n)-diag(ones(1,n)));
A_xi = (1-xi)*A_S + xi*A_C;

end
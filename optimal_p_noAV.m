function cvx_optval = optimal_p_noAV(n,xi,beta)
A = S2C_matrix(n,xi);
% A = [0,0.2,0.8;0.4,0,0.6;0,1,0];

cvx_begin quiet
    variables p(n,1) delta(n) y(n,n)
    maximize(sum(p.*(1-p))-sum(delta))
    subject to 
        sum(y,2) == transpose(beta*(sum(A.*((1-p)*ones(1,n)),1)+sum(y,1)))+delta-(1-p)
        p >= 0
        delta >= 0
        y >= 0
        p <= 1
cvx_end    
end
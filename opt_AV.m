function solution = opt_AV(n,xi,beta,k)
A = S2C_matrix(n,xi);
cvx_begin quiet
    variables p(n,1) delta(n,1) x(n,1) r(n,n) z(n,1)
    maximize( sum(p.*(1-p),1)-sum(delta,1)-k*sum(z,1) )
    subject to 
        x == beta*(transpose(sum(A.*(x*ones(1,n)))))+delta
        z == transpose(sum(A.*((1-p-x)*ones(1,n)),1)+sum(r,1))
        sum(r,2) == z - (1-p-x)   
%         p>=0
        delta>=0
        z>=0
        r>=0
        x>=0
cvx_end
solution.profit = cvx_optval;
solution.price = p;
solution.delta = delta;
solution.z = z;
solution.r = r;
solution.x = x;

end
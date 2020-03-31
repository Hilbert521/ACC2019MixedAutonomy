function solution = AV_gu(n,xi,beta,k)
A = S2C_matrix(n,xi);

% A = [0,0.2,0.8;0.4,0,0.6;0,1,0];

p= sdpvar(n,1) ;
delta= sdpvar(n,1); 
y = sdpvar(n,n,'full');
r = sdpvar(n,n,'full');
x = sdpvar(n,1);
z = sdpvar(n,1);
d = sdpvar(n,1);

constraints=[d==min(x+z,(1-p)),...
    x'==beta*((min(x,1-p))'*A+ones(1,n)*y)+delta',...
    y*ones(n,1) == max((x-(1-p)),0),...
    z'==min(z,max(1-p-x,0))'*A+ones(1,n)*r,...
    r*ones(n,1)==max(z-max(((1-p)-x),0),0),...
    r(1,2) == r(1,3),...
    y(1,2) == y(1,3),...
    delta(2) == delta(3),...
    delta>=0,...
    y(:)==0,...
    r(:)>=0,...   
    p>=0,...
    x==0,...
    z>=0,...
    x+z>=(1-p),... %This needs to hold for the assumption in the cost to be correct
    ];

obj=sum(p.*(1-p))-sum(delta)-k*sum(z);
options = sdpsettings('verbose',0,'solver','gurobi');

sol = optimize(constraints,-obj,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution.profit = value(obj);
 solution.price = value(p);
 solution.delta = value(delta);
 solution.z = value(z);
 solution.r = value(r);
 solution.y = value(y);
 solution.x = value(x);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem);
end

end
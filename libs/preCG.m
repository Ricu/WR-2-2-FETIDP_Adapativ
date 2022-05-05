function [x,resid,iter,kappa_est,alpha,beta] = preCG(A,L,b,x0,tol,ploth)
invM = @(x) L'\(L\x);
rk = b - A(x0);
zk = invM(rk);
pk = zk;
xk = x0;
iter = 0;
alpha_vec = zeros(1000,1);
beta_vec = zeros(1000,1);
normb = norm(b);


while norm(rk)/normb > tol
    if nargin > 5 && iter < 3
        ploth(xk,iter);
    end
    ak = (rk'*zk) / (pk'*A(pk));
    xk = xk + ak * pk;
    rkp1 = rk - ak * A(pk);
    zkp1 = invM(rkp1);
    bk = (zkp1' * rkp1)/(zk' * rk);
    pk = zkp1 + bk * pk;
    rk = rkp1;
    zk = zkp1;
    iter = iter+1;
    if iter > size(alpha_vec)
        alpha_vec = [alpha_vec ; zeros(1000,1)];
        beta_vec = [beta_vec ; zeros(1000,1)];
    end
    alpha_vec(iter) = ak;
    beta_vec(iter) = bk;
    
%     fprintf("R= %e, RR= %e, iter=%i\n",norm(rk),norm(rk)/normb,iter);
end

x = xk;
resid = rk;
alpha = alpha_vec(1:iter);
beta = beta_vec(1:iter);

if nargout > 3
    temp1 = [sqrt(beta(1:end-1))./alpha(1:end-1); 0];
    temp2 = (1./alpha) + [0;beta(1:iter-1)]./[1;alpha(1:iter-1)];
    temp3 = [0; sqrt(beta(1:end-1))./alpha(1:end-1)];
    Tk = spdiags([temp1 temp2 temp3],-1:1,iter,iter);
    kappa_est = cond(full(Tk)); 
end
end


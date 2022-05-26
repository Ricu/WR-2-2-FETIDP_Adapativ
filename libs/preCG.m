function [x,resid,iter,kappa_est,alpha,beta] = preCG(A,invM,b,x0,tol,VK,ploth,U,invUFU,d)

r0 = b - A(x0);
rk = r0;
z0 = invM(rk);
zk = z0;
pk = zk;
xk = x0;
iter = 0;
alpha_vec = zeros(1000,1);
beta_vec = zeros(1000,1);
nanInfo = cell(1000,1);

figure("Name","Loesungen waehrend der Iteration von PCG")
while norm(zk)/norm(z0) > tol %norm(rk)/norm(r0) > tol %
    if nargin > 5 && iter < 4
        if strcmp('Deflation',VK)
            xBar = U*invUFU*U'*d;   % Korrektur bei Deflation-VK notwendig
            ploth(xk+xBar,iter,VK);
        else
            ploth(xk,iter,VK);
        end
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
    
    nanInfo{iter} = isnan(xk);
%     fprintf('Iter %i, Residuum %10f, Anzahl NaN %i\n',iter,norm(rk)/norm(r0),nnz(isnan(xk)))
end

x = xk;
resid = rk;
alpha = alpha_vec(1:iter);
beta = beta_vec(1:iter);
nanInfo = nanInfo(1:iter);

if nargout > 3
    temp1 = [sqrt(beta(1:end-1))./alpha(1:end-1); 0];
    temp2 = (1./alpha) + [0;beta(1:iter-1)]./[1;alpha(1:iter-1)];
    temp3 = [0; sqrt(beta(1:end-1))./alpha(1:end-1)];
    Tk = spdiags([temp1 temp2 temp3],-1:1,iter,iter);
    if ~strcmp('Deflation',VK)
        kappa_est = cond(full(Tk));
    else
        kappa_est = 0;
    end
end
end


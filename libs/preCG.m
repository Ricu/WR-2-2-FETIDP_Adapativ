function [x,iter,kappa_est,residual] = preCG(A,invM,b,pcg_param,VK,plot_struct,constraint_struct)
% Input: A: Systemmatrix
% Input: invM: Vorkonditionierer (function handle)
% Input: b: Rechte Seite 
% Input: pcg_param: Structure mit PCG-Parametern
%        Komponenten: tol, x0 (function handle), resid_type
% Input: VK: String mit Name des VK
% Input: plot_struct: Structure mit Plot-Komponenten
%        Komponenten: plot_iteration(boolean),ploth(function handle)
% Input: constraint_struct: Structure mit allen Parametern fuer constraints 
%        Komponenten: U,invUFU,IminusPtranspose


% Output: x: Loesungsvektor des LGS
% Output: iter: Anzahl Iterationen
% Output: kappa_est: Konditionszahlschaetzung aus dem Lanczos-Prozess
% Output: residual: Vektor mit Residuen aus der Iteration

%% Entpacken der Structures
tol          = pcg_param.tol;
x0           = pcg_param.x0(length(b));
resid_type   = pcg_param.resid_type;

% Definiere Korrektur Matrix (fuer Loesung mit Deflation)
correction_matrix = constraint_struct.U*constraint_struct.invUFU*constraint_struct.U';
IminusPtranspose = constraint_struct.IminusPtranspose;

%% Initialisierung
xk = x0;        % Loesungsvektor
r0 = b - A(x0); % Residuum
rk = r0;
z0 = invM(rk);  % Vorkonditioniertes Residuum
zk = z0;
pk = zk;        % Abstiegsrichtung

%Speicherreservierung fuer alpha, beta und die Residuen
alpha_vec = zeros(1000,1);
beta_vec = zeros(1000,1);
residual_vec = zeros(1000,1);

% Definiere betrachtetes Residuum fuer Abbruchbedingung
if strcmp('vorkonditioniert',resid_type)
    residual = norm(zk)/norm(z0);
elseif strcmp('Deflation',VK) && strcmp('nicht-vorkonditioniert,alternativ',resid_type)
    residual = norm(IminusPtranspose(rk))/norm(IminusPtranspose(r0));    
else % nicht-vorkonditioniert
    residual = norm(rk)/norm(r0);
end

if plot_struct.plot_iteration
    figure("Name","Loesungen waehrend der Iteration von PCG")
end

%% Iteration bis geforderte Genauigkeit erreicht
iter = 0;   % Anzahl Iterationen
while residual > tol     % Abbruchbedingung
    % Plotten der Loesungen der ersten Iterationen
    if plot_struct.plot_iteration && iter < 4
        if strcmp('Deflation',VK) && iter > 0
            % Korrektur der Loesung bei Deflation-VK notwendig
            xBar = correction_matrix*b;
            plot_struct.ploth(xk+xBar,iter,VK);
        else
            plot_struct.ploth(xk,iter,VK);
        end
    end
    
    % Verfahren
    ak = (rk'*zk) / (pk'*A(pk));    
    xk = xk + ak * pk;              % Loesungsvektor
    rkp1 = rk - ak * A(pk);         % Neues Residuum
    zkp1 = invM(rkp1);              % Vorkonditioniertes neues Residuum
    bk = (zkp1' * rkp1)/(zk' * rk);
    pk = zkp1 + bk * pk;            % Abstiegsrichtung
    
    % Aktualisiere Parameter fuer naechste Iteration
    rk = rkp1;
    zk = zkp1;
    iter = iter+1;
    if strcmp('vorkonditioniert',resid_type)
        residual = norm(zk)/norm(z0);
    elseif strcmp('Deflation',VK) && strcmp('nicht-vorkonditioniert,alternativ',resid_type)
        residual = norm(IminusPtranspose(rk))/norm(IminusPtranspose(r0));
    else % nicht-vorkonditioniert
        residual = norm(rk)/norm(r0);
    end
    
    % Speichern der berechneten Werte
    if iter > size(alpha_vec)
        alpha_vec = [alpha_vec ; zeros(1000,1)];
        beta_vec = [beta_vec ; zeros(1000,1)];
        residual_vec = [residual_vec ; zeros(1000,1)];
    end    
    alpha_vec(iter) = ak;
    beta_vec(iter) = bk;
    residual_vec(iter) = residual;
end

% Aufstellen des Loesungsvektors
x = xk;

% Entferne nicht benoetigten Speicherplatz
alpha = alpha_vec(1:iter);
beta = beta_vec(1:iter);
residual = residual_vec(1:iter);

% Berechne die Konditionszahlschätzung des Lanczos-Prozess
if nargout > 3
    temp1 = [sqrt(beta(1:end-1))./alpha(1:end-1); 0];
    temp2 = (1./alpha) + [0;beta(1:iter-1)]./[1;alpha(1:iter-1)];
    temp3 = [0; sqrt(beta(1:end-1))./alpha(1:end-1)];
    Tk = spdiags([temp1 temp2 temp3],-1:1,iter,iter);
    if nnz(isnan(x)) > 0
        kappa_est = 0;
        fprintf('Die Konditionszahl beim %s-Vorkonditionierer konnte nicht berechnet werden, da die Loesung NaN enthaelt.\n',VK)
    else
        kappa_est = cond(full(Tk));
    end
end
end




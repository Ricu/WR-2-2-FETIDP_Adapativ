function [phi,d_phi] = baseFun(order)
if order == 1
    %% P1
        phi1 = @(x,y) 1-x-y; phi2 = @(x,y) x; phi3 = @(x,y) y; 
        phi = {phi1,phi2,phi3};
        d_phi = {@(x,y) -1,@(x,y) 1, @(x,y) 0;  %dx
               @(x,y) -1,@(x,y) 0, @(x,y) 1};   %dy
elseif order == 2
    %% P2
    vertices = [0,0;
                1,0;
                0,1;
                1/2,0;
                1/2,1/2;
                0,1/2];
    numBaseFun = 6;
    phi = cell(1,numBaseFun);
    % Zeile 1: dx   Zeile 2:  dx
    d_phi = cell(2,numBaseFun);
    A = zeros(numBaseFun);
    A(:,1) = ones(numBaseFun,1);
    A(:,2) = vertices(:,1);
    A(:,3) = vertices(:,2);
    A(:,4) = vertices(:,1).^2;
    A(:,5) = vertices(:,1).*vertices(:,2);
    A(:,6) = vertices(:,2).^2;
    for i = 1:numBaseFun
        b = zeros(numBaseFun,1);
        b(i) = 1;

        c = A\b;
        phi{i} = @(x,y) c(1) + c(2)*x + c(3)*y + c(4)*x.^2 +c(5)*x.*y + c(6)*y.^2;
        %dx
        d_phi{1,i}  = @(x,y) +c(2) + 2*c(4)*x + c(5)*y;
        %dy
        d_phi{2,i}  = @(x,y) +c(3)+ c(5)*x + 2*c(6)*y;
    end
end
end
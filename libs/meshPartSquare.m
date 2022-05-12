function [x__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,x,tri)
% Gitterpartitionierung von [0,1]^2 in N x N Teilquadrate.
%
% INPUT:
%    N^2: Anzahl Teilgebiete
%    x:   Punkteliste
%    tri: Elementliste (Dreiecke)
%
% OUTPUT:
%    x__sd:   Cell array; lokale Punkteliste
%    tri__sd: Cell array; lokale Elementliste
%    l2g__sd: Cell array; local-to-global-map
%    logicalTri__sd: Cell array; logische Liste, welche Dreiecke in welchem
%    TG sind
   %
    % Berechne den Schwerpunkt der Dreiecke.
    c = (x(tri(:,1),:) + x(tri(:,2),:) + x(tri(:,3),:))/3;
    
    % Aufteilung von [0,1] mit aequidistanten Stuetzstellen.
    lin = linspace(0,1,N+1);
    
    x__sd = cell(N^2,1);
    tri__sd = cell(N^2,1);
    l2g__sd = cell(N^2,1);
    logicalTri__sd  = cell(N^2,1);
    cnt = 0;
    
    for i = 1:N
        for j = 1:N
            cnt = cnt + 1;
            
            % Finde alle Dreiecke im Teilgebiet (i,j).
            b = (c(:,1) >= lin(j)-1e-13) & (c(:,1) <= lin(j+1)+1e-13) & ...
                (c(:,2) >= lin(i)-1e-13) & (c(:,2) <= lin(i+1)+1e-13);
            logicalTri__sd{cnt} = b;    
            
            tri__sd{cnt} = tri(b,:); % globale Nummerierung
            
            % Finde alle Knoten im Teilgebiet (i,j).
            b = false(size(x,1),1);
            b(tri__sd{cnt}(:)) = true;
            x__sd{cnt} = x(b,:);
            l2g__sd{cnt} = find(b);
            
            % Mappe globale Nummerierung der lokalen Triangulierung 
            % auf lokale Nummerierung.
            map = zeros(size(x,1),1);
            map(b) = 1:nnz(b);
            tri__sd{cnt} = map(tri__sd{cnt});
        end
    end
end
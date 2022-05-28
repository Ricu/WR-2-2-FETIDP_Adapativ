function [rhoTri,rhoTriSD,indElementsrhoMax,maxRhoVert,maxRhoVertSD] = coefficient_2iii(rhoMax_vec,rhoMin,tri,numVert,numTri,numSD,logicalTri__sd)

%% Definiere Koeffizientenfunktion auf den Elementen
% Setze alle Koeffizienten der Elemente auf 'rhoMin'
rhoTri = rhoMin*ones(numTri,1);
for r = 1:length(rhoMax_vec)
    rhoMax = rhoMax_vec(r);
    % Setze nun alle 'rhoMax' Koeffizienten.
    rhoTri(rand(length(rhoTri),1) < 0.25) = rhoMax;
end
indElementsrhoMax = (rhoTri == rhoMax); % Logischer Vektor, welche Elemente in rhoMax liegen

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); 
maxRhoVertSD = cell(numVert,1);
for i = 1:numVert % Iteriere ueber Knoten
    iVec = i*ones(1,size(tri,2));
    cnt = 1;
    for j = 1:numTri % Iteriere ueber Dreiecke
        testMembership = ismember(iVec,tri(j,:)); 
        if nnz(testMembership) > 1    % Pruefe, ob Knoten im Dreieck liegt
            vertTris{i}(cnt) = j;   % Enthaelt fuer jeden Knoten die Dreiecke in denen er liegt
            cnt = cnt+1;
        end
    end
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten teilgebietsweise
    for k = 1:numSD
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % logischer Vektor welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % maximalen Koeffizienten pro Knoten teilgebietsweise
    end
end
end

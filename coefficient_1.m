function [rhoTri,rhoTriSD,indElementsCanal,maxRhoVert,maxRhoVertSD] = coefficient(xMin,xMax,yMin,yMax,rhoCanal,rhoNotCanal,vert,tri,numVert,numTri,numSD,logicalTri__sd)

%% Definiere Koeffizientenfunktion auf den Elementen
indVertCanal = (xMin <= vert(:,1) & vert(:,1) <= xMax & yMin <= vert(:,2) & vert(:,2) <= yMax);  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = find(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen

rhoTri = rhoNotCanal*ones(numTri,1); % Außerhalb des Kanals entspricht die Koeffizientenfunktion 1
for i=1:numTri % Iteriere ueber die Elemente
    if ismember(tri(i,:),numVertCanal) % Alle Knoten des Elements liegen im Kanal
        rhoTri(i)=rhoCanal;    % Im Kanal entspricht die Koeffizientenfunktion 10^6
    end
end
indElementsCanal = rhoTri > 1; % Logischer Vektor, welche Elemente im Kanal liegen

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


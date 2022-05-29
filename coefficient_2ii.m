function [rhoTri,rhoTriSD,indElementsrhoMax,maxRhoVert,maxRhoVertSD] = coefficient_2ii(rhoMax_vec,rhoMin,rhoMax_sd,l2g__sd,tri,numVert,numTri,numSD,logicalTri__sd)
% Input: 
% Input: 
% Input: tri: Elementliste
% Input: numVert,numTri,numSD: Anzahl Knoten, Elemente, Teilgebiete
% Input: logicalTri__sd: Logischer Vektor, welche Dreiecke in welchem TG enthalten sind

% Output: rhoTri,rhoTriSD: Koeffizient pro Element (und teilgebietsweise)
% Output: indElementsrhoMax: Logischer Vektor, welche Elemente in rhoMax liegen
% Output: maxRhoVert,maxRhoVertSD: maximaler Koeffizient pro Knoten (und teilgebietsweise)

%% Definiere Koeffizientenfunktion auf den Elementen
rhoTri = rhoMin*ones(numTri,1); % Die Koeffizientenfunktion entspricht rhoMin außerhalb der markierten TG

numVertrhoMax = []; % Knotennummern der Knoten, die in rhoMax liegen
for i = 1: length(rhoMax_sd)
    numVertrhoMax = [numVertrhoMax l2g__sd{rhoMax_sd(i),1}];
end
for r = 1:length(rhoMax_vec)
    rhoMax = rhoMax_vec(r);
    for i=1:numTri % Iteriere ueber die Elemente
        if ismember(tri(i,:),numVertrhoMax) % Alle Knoten des Elements liegen im Kanal
            rhoTri(i)=rhoMax;    % Im Kanal entspricht die Koeffizientenfunktion 10^6
        end
    end
end
indElementsrhoMax = (rhoTri == rhoMax); % Logischer Vektor, welche Elemente in rhoMax liegen

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i}); % Koeffizientenfunktion pro Element teilgebietsweise
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); 
maxRhoVertSD = cell(numVert,1);
for i = 1:numVert % Iteriere ueber Knoten
    cnt = 1;
    for j = 1:numTri % Iteriere ueber Dreiecke  
        if ismember(i,tri(j,:)) % Pruefe, ob Knoten im Dreieck liegt
            vertTris{i}(cnt) = j;   % Enthaelt fuer jeden Knoten die Dreiecke in denen er liegt
            cnt = cnt+1;
        end
    end
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten teilgebietsweise
    for k = 1:numSD % Iteriere ueber TG
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % Logischer Vektor, welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % Maximaler Koeffizient pro Knoten teilgebietsweise
    end
end
end


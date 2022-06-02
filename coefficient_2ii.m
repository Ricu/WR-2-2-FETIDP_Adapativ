function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_2ii(rhoMax,rhoMin,rhoMax_sd,l2g__sd,tri,vert,numVert,numTri,numSD,logicalTri__sd,N,plot)
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

for ind = 1:length(rhoMax_sd)
    sd = rhoMax_sd(ind);
    rhoTri(logicalTri__sd{sd}) = rhoMax;
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
    [vertTris{i},~,~] = find(i == tri);
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten teilgebietsweise
    for k = 1:numSD % Iteriere ueber TG
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % Logischer Vektor, welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % Maximaler Koeffizient pro Knoten teilgebietsweise
    end
end

if plot == true
    %% Plotten des Gitters mit Kanal
    figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
    patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]); hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(indElementsrhoMax,:),'edgecol','k','facecol',[.8,.9,1]);
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1, 'color', 'r')
    end
    legend('\rho = 1','\rho = 10^6','Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
end

end


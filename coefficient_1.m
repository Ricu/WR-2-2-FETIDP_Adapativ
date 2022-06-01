function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_1(xMin,xMax,yMin,yMax,rhoCanal,rhoNotCanal,vert,tri,numVert,numTri,numSD,logicalTri__sd,N,plot)
% Input: xMin,xMax,yMin,yMax: Grenzen des Kanalgebiets in x- und y-Richtung
% Input: rhoCanal,rhoNotCanal: rho im Kanal und außerhalb des Kanals
% Input: vert,tri: Knoten- und Elementliste
% Input: numVert,numTri,numSD: Anzahl Knoten, Elemente, Teilgebiete
% Input: logicalTri__sd: Logischer Vektor, welche Dreiecke in welchem TG enthalten sind

% Output: rhoTri,rhoTriSD: Koeffizient pro Element (und teilgebietsweise)
% Output: indElementsCanal: Logischer Vektor, welche Elemente im Kanal liegen
% Output: maxRhoVert,maxRhoVertSD: maximaler Koeffizient pro Knoten (und teilgebietsweise)

%% Definiere Koeffizientenfunktion auf den Elementen
indVertCanal = (xMin <= vert(:,1) & vert(:,1) <= xMax & yMin <= vert(:,2) & vert(:,2) <= yMax);  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = find(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen

rhoTri = rhoNotCanal*ones(numTri,1); % Koeffizient (zuerst) auf allen Elementen = Koeffizient (gleich) auf Elementen außerhalb des Kanals
for i=1:numTri % Iteriere ueber die Elemente
    if ismember(tri(i,:),numVertCanal) % Alle Knoten des Elements liegen im Kanal und damit das Element selber
        rhoTri(i)=rhoCanal;    % Koeffizient auf Elementen innerhalb des Kanals
    end
end
indElementsCanal = rhoTri > 1; % Logischer Vektor, welche Elemente im Kanal liegen

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});  % Koeffizientenfunktion pro Element teilgebietsweise
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
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % Maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten eines TG
    for k = 1:numSD % Iteriere ueber TG
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % Logischer Vektor, welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % Maximaler Koeffizient pro Knoten teilgebietsweise
    end
end

if plot == true
    %% Plotten des Gitters mit Kanal
    figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
    patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1],'edgecolor',"#3c3c3c"); 
    hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(indElementsCanal,:),'edgecol','k','facecol',"#2b8cbe",'edgecolor',"#3c3c3c");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
    rhoMax = sprintf('\\rho = %i',rhoCanal);
    legend('\rho = 1',rhoMax,'Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
end

end


function cEdges = findEdges(cDual,dual,l2g__sd)
% Finde Kanten (alle Eckknoten sind primal).
%
% INPUT:
%    cDual: Cell-array: Liste (fuer jedes Teilgebiet) aller lokaler dualer 
%           Knoten mit lokaler Nummerierung (cDual{i} darf ein logischer 
%           Vektor sein oder eine Nummernliste).
%    dual: Logischer Vektor: Globale Liste aller dualer Knoten.
%    l2g__sd: Cell-array: Local-2-global-map fuer KnotenIDs.
%
% OUTPUT:
%    cEdges: Cell-array: Liste aller Kanten; globale Knotennummerierung.
%
    %
    
    % Bestimme Teilgebietszugehoerigkeit dualer Knoten / Kantenknoten.
    dual_sd = findAdjacentSubdomainsOfDualNodes(cDual,l2g__sd,dual);
    
    % Sortiere Teilgebietsnummern bzgl. der Spalten, sodass
    %    dual_sd(:,1) < dual_sd(:,2).
    dual_sd = sort(dual_sd,2);
    
    % Sortiere Zeilen.
    % [ 1,2   <-- Kante 1
    %   1,2   <-- Kante 1
    %   1,2   <-- Kante 1
    %   1,3   <-- Kante 2
    %   1,3   <-- Kante 2
    %   ...
    %   3,4   <-- Kante m
    %   3,4]  <-- Kante m
    [dual_sd,sortmap] = sortrows(dual_sd);
    
    % Ziehe Zeilen voneinander ab (von oben nach unten).
    % [ 1,2   <-- Kante 1  --|
    %   1,2   <-- Kante 1  <-|  --|
    %   1,2   <-- Kante 1       <-|  --|
    %   1,3   <-- Kante 2  --|       <-|
    %   1,3   <-- Kante 2  <-|  --|
    %   ...                     <-|  --|
    %   3,4   <-- Kante m  --|       <-|
    %   3,4]  <-- Kante m  <-|
    % Ergibt:
    % [ 0,0
    %   0,0
    %   0,1
    %   0,0
    %   ...
    %   ?,?
    %   0,0]
    tmp = diff(dual_sd);
    
    % Alle Zeilen, die vollstaendig Null sind, werden auf 'false' 
    % gesetzt (true --> eindeutige Zeile).
    % [ 0,0  --> false
    %   0,0  --> false
    %   0,1  --> true
    %   0,0  --> false
    %   ...
    %   ?,?  --> ?
    %   0,0]  --> false
    tmp = any(tmp,2);
    
    % Die erste Zeile von 'dual_sd' ist per Definition eindeutig.
    % Fuege also 'true' vorne ein.
    ind = [true;
           tmp];
    
    % Fuege am Ende noch ein dummy-'true' ein; eine Art Stoppmarker.
    ind = [ind;
           true];
    
    % Suche nun alle eindeutigen Zeilen.
    % find(indS) = [1;3;...]
    % dual_sd(indS(1:end-1),:) = [1,2
    %                             1,3
    %                             ...
    %                             3,4]
    ind_unique = find(ind);
    nEdges = length(ind_unique)-1;
    
    % Bestimme globale KnotenIDs der dualen Knoten und sortiere diese 
    % genauso wie 'dual_sd'.
    nodesDual = find(dual);
    nodesDual = nodesDual(sortmap);
    
    % Wir koennen nun die Kantenknoten extrahieren.
    cEdges = cell(nEdges,1);
    for i = 1:nEdges
        cEdges{i} = nodesDual(ind_unique(i):ind_unique(i+1)-1);
    end
end

function dual_sd = findAdjacentSubdomainsOfDualNodes(cDual,l2g__sd,dual)
    % Bestimme Teilgebietszugehoerigkeit.
    n = length(dual);
    column_cnt = zeros(n,1);
    dual_sd = zeros(n,2);  % zwei Teilgebiete pro dualem Knoten
    for i = 1:length(l2g__sd)
        nodeIDs = l2g__sd{i}(cDual{i});  % globale KnotenIDs dualer Knoten
        column = column_cnt(nodeIDs)+1;
        for j = 1:length(column)
            dual_sd(nodeIDs(j),column(j)) = i;
        end
        column_cnt(nodeIDs) = column_cnt(nodeIDs) + 1;
    end
    dual_sd = dual_sd(dual,:);
end
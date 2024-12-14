function [G,B,I] = compute_GandB(h,d,A,reservoirsID,E,pipes_in_zone)

N=length(h);
% E = d.NodesConnectingLinksIndex;
nEdges = size(E,1);

% Roughness (unitless) %

roughness = d.LinkRoughnessCoeff(pipes_in_zone)';
% roughness(end-3:end) = 140;

% Diameters (mm in EPANET) %

diameters = d.LinkDiameter(pipes_in_zone)'/1000; % Diameter (mm to m)

% Diameters (m in EPANET) %

lengths = d.LinkLength(pipes_in_zone)'; % Length (d.LinkLength in m)
% lengths(end-3:end) = 1;

G=zeros(N);%((roughness.^1.852).*(diameters.^4.8704))./(10.675*lengths);

for i=1:nEdges
    G(E(i,1),E(i,2)) = ((roughness(i)^1.852)*(diameters(i)^4.8704))/(10.675*lengths(i));
    G(E(i,2),E(i,1)) = G(E(i,1),E(i,2));
end

if nargin > 1

    B = zeros(N); I = zeros(nEdges,N);
    for i=1:length(E)
        if h(E(i,1)) >= h(E(i,2))
            B(E(i,1),E(i,2)) = 1; B(E(i,2),E(i,1)) = -1;
            I(i,E(i,1)) = 1; I(i,E(i,2)) = -1;
        else
            B(E(i,1),E(i,2)) = -1; B(E(i,2),E(i,1)) = 1;
            I(i,E(i,1)) = -1; I(i,E(i,2)) = 1;
        end
    end

else

    % Aproximated incidence matrix computation %

    innernodesIDs = setdiff(1:N,reservoirsID);
    C = zeros(size(A));
    for i=1:length(reservoirsID)
        for j =1:length(innernodesIDs)
            path = shortestpath(G,reservoirsID(i),innernodesIDs(j));
            for k=1:length(path)-1
                C(path(k),path(k+1)) = C(path(k),path(k+1)) + 1;
            end
        end
    end

    B = zeros(nEdges,N);
    for i=1:length(E)
        if C(E(i,1),E(i,2)) > C(E(i,2),E(i,1))
            B(E(i,1),E(i,2)) = -1; B(E(i,2),E(i,1)) = 1;
        else
            B(E(i,1),E(i,2)) = 1; B(E(i,2),E(i,1)) = -1;
        end
    end

    B = -B;

end

end
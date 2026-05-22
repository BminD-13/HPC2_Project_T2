function [u,energy] = solveLame(coordinates, elements, dirichlet, neumann, f, g, uD, lambda, mu)

nC = size(coordinates,1);

% globale Matrix + RHS
A = sparse(2*nC, 2*nC);
b = zeros(2*nC,1);

%*** Assembly (nur Dreiecke)
for j = 1:size(elements,1)
    
    nodes = elements(j,:);
    coordsT = coordinates(nodes,:);
    
    % DOF-Indizes (2 pro Knoten!)
    I = 2*nodes([1 1 2 2 3 3]) - [1 0 1 0 1 0];
    
    % Elementmatrix (brauchst du!)
    A(I,I) = A(I,I) + stima3(coordsT, lambda, mu);
    
    % Volumenkraft
    xm = mean(coordsT,1);
    fs = f(xm)';   % -> Spaltenvektor
    
    area = det([1 1 1; coordsT'])/2;
    
    b(I) = b(I) + area * [fs; fs; fs] / 3;
end

%*** Neumann
if ~isempty(neumann)
    for j = 1:size(neumann,1)
        nodes = neumann(j,:);
        coordsE = coordinates(nodes,:);
        
        I = 2*nodes([1 1 2 2]) - [1 0 1 0];
        
        xm = mean(coordsE,1);
        gm = g(xm)'; 
        
        edgeLength = norm(coordsE(2,:) - coordsE(1,:));
        
        b(I) = b(I) + edgeLength * [gm; gm] / 2;
    end
end

%*** Dirichlet
DirichletNodes = unique(dirichlet);
[W,M] = uD(coordinates(DirichletNodes,:));

B = sparse(size(W,1),2*nC);

for k = 0:1
    for l = 0:1
        B(1+l:2:size(M,1),2*DirichletNodes-1+k) = ...
            diag(M(1+l:2:size(M,1),1+k));
    end
end

mask = find(sum(abs(B),2));
A = [A, B(mask,:)'; B(mask,:), sparse(length(mask),length(mask))];
b = [b; W(mask)];

%*** lösen
x = A \ b;
u = x(1:2*nC);

% Energie (optional)
energy = u' * A(1:2*nC,1:2*nC) * u;

end
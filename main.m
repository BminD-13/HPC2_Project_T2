clc
clear 
close all

%% Geometrie

% Traktor vernetzt 
mesh

% Kreis
theta = linspace(0,2*pi,200);
circle = [3 + 2*cos(theta); 3 + 2*sin(theta)];

%% Randbedingungen

f  = @(x) zeros(size(x,1),1);
g  = @(x) zeros(size(x,1),1);

% Segment-Knoten
nodes5 = unique(boundary{6}(:));
nodes3 = unique(boundary{13}(:));
nodes7 = unique([boundary{14}(:);boundary{15}(:)]);

% Zugehörige Koordinaten
coords5 = coordinates(nodes5,:);
coords3 = coordinates(nodes3,:);
coords7 = coordinates(nodes7,:);

uD = @(x) uDiriclet(x, coords5, coords3, coords7); % am Ende des Skripts

%% PDE Solver

[x,Energie] = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD);

%% Integral über den Kreis



% Chat GPT Anwort:

% Bei P1-FEM (lineare Dreieckselemente) gilt auf jedem Dreieck:
% lineare Basisfunktionen (hat functions, baryzentrische Koordinaten)

% Das nennt man:
% Finite-Element-Rekonstruktion
% Finite-Element-Interpolation
% P1-Interpolant
% Piecewise linear interpolation (stückweise lineare Interpolation)

%% Plot

figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), x);

axis equal;
view(3);
shading interp;
colorbar;

title('P1-FEM Lösung der Laplace-Gleichung');
xlabel('x');
ylabel('y');
zlabel('u_h');

disp(['Energie: ', num2str(Energie)]);

%% Funktionen

function W = uDiriclet(x,coords5,coords3,coords7)

tol = 1e-10;

n = size(x,1);

W = zeros(n,1);

for i = 1:n

    if any(vecnorm(coords5 - x(i,:),2,2) < tol)

        W(i) = 5;

    elseif any(vecnorm(coords3 - x(i,:),2,2) < tol)

        W(i) = 3;

    elseif any(vecnorm(coords7 - x(i,:),2,2) < tol)

        W(i) = 7;

    end
end

end
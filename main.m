clc
clear 
close all

%% Geometrie

% Traktor vernetzt 
mesh

% Kreis
radius = 2;
cx = 2.5; cy = -1.5;
N = 2000;
theta = linspace(0,2*pi,N);
circle = [cx + radius*cos(theta); cy + radius*sin(theta)]';

%% Randbedingungen

f  = @(x) zeros(size(x,1),1);
g  = @(x) zeros(size(x,1),1);

% Segment-Knoten
nodes5 = unique(boundary{6}(:));
nodes3 = unique(boundary{13}(:));
nodes7 = unique([boundary{14}(:);boundary{15}(:)]);
coords5 = coordinates(nodes5,:);
coords3 = coordinates(nodes3,:);
coords7 = coordinates(nodes7,:);

uD = @(x) uDiriclet(x, coords5, coords3, coords7);

%% PDE Solver

[x,Energie] = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD);

%% Integral über den Kreis

ds = 2*pi*radius / (N-1);

summe = 0;

for i = 1:N-1
    p = circle(i,:);
    [P,U] = findTriangle(p, coordinates, elements, x);
    update = evalP1(p(1),p(2),P,U);
    summe = summe + update;
end

integral = summe * ds; % 54-81

% Analytische Abschätzung:
% MinMax Werte laut Grafik sind ca. 4.3 und 6.5 und Radius = 2
% Imin = 2 * r * pi * min = 2 * 2 * pi * 4.3 = 54
% Imax = 2 * r * pi * max = 2 * 2 * pi * 6.5 = 81


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

% Kreis
hold on
plot(circle(:,1), circle(:,2), 'r', 'LineWidth', 2)
hold off


disp(['Energie: ', num2str(Energie)]);

%% Funktionen

% Randbedingung
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

% P1 Interpolation

function inside = pointInTriangle(p, P)

A = [
    1      1      1;
    P(1,1) P(2,1) P(3,1);
    P(1,2) P(2,2) P(3,2)
];

rhs = [1; p(1); p(2)];

lambda = A \ rhs;

inside = all(lambda >= -1e-12) && all(lambda <= 1+1e-12);

end


function val = evalP1(x,y,P,U)

A = [
    1 1 1;
    P(:,1)';
    P(:,2)'
];

rhs = [1; x; y];

lambda = A \ rhs;

val = lambda' * U;

end

function [P,U] = findTriangle(p, coordinates, elements, u)

for k = 1:size(elements,1)

    nodes = elements(k,:);
    Pk = coordinates(nodes,:);
    
    if pointInTriangle(p, Pk)

        P = Pk;
        U = u(nodes);
        return;

    end
end

error('Point not in mesh domain');
end

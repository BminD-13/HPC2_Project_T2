clc
clear 
close all

[coords, elements, boundary] = createGeometry(0.5);

dirichlet = boundary{1}
neumann = boundary{2}

f = @(x) ones(size(x,1),1)*0.01;
g = @(x) zeros(size(x,1),1);
uD = @(x) zeros(size(x,1),1);

[x,Energie] = solveLaplace(coords, elements, dirichlet, neumann, f, g, uD);

figure;
trisurf(elements, coords(:,1), coords(:,2), x);

axis equal;
view(3);
shading interp;
colorbar;

title('P1-FEM Lösung der Laplace-Gleichung');
xlabel('x');
ylabel('y');
zlabel('u_h');

disp(['Energie: ', num2str(Energie)]);

function [int_u] = integral_u(circle, radius, x, coordinates, elements)

N = length(circle);
ds = 2*pi*radius / (N-1);
u_along_circle = [];

for i = 1:N-1
    p       = circle(i,:);
    [P,U]   = findTriangle(p, coordinates, elements, x);
    plane_coefficients = planeP1(P,U);
    point   = [1,p]';
    
    % u
    u = plane_coefficients' * point; % [a b c] * [1 x y]'
    u_along_circle(end+1) = u;
end

int_u = sum(u_along_circle) * ds; % in [54,81]
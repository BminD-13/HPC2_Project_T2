function [int_dudn] = integral_dudn(circle, center, radius, x, coordinates, elements)
N = length(circle);
ds = 2*pi*radius / (N-1);
dudn_along_circle = [];

for i = 1:N-1
    p       = circle(i,:);
    [P,U]   = findTriangle(p, coordinates, elements, x);
    plane_coefficients = planeP1(P,U);

    normal = (p-center) / norm(p-center, 2);
    dudn = plane_coefficients(2:3)' * normal'; % [b c] * [mx,my]'
    dudn_along_circle(end+1) = dudn;
end

int_dudn = sum(dudn_along_circle) * ds;

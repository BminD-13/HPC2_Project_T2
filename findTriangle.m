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
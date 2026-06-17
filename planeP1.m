function coeff = planeP1(P,U)
    M = [
        ones(3,1), P
    ];
    coeff = M \ U;
end


function [coordinates,elements,boundary] = createGeometry(hmax)

    g = [
     % AB  BC  CM  MD  DE  EF  FG  GH  HI  IJ  JK  KL  LA  MM  MM  NN  NN;
        1,  2,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  1,  1,  1;
        2,  8, 10, 14, 16, 16, 14, 14, 13, 13,  8,  8,  4,  4,  6, 11, 13;
        8, 10, 14, 16, 16, 14, 14, 13, 13,  8,  8,  4,  2,  6,  4, 13, 11; 
        2,  2,  2,  2,  2, 10, 10, 16, 16, 10, 14, 20, 20,  2,  2,  2,  2;
        2,  2,  2,  2, 10, 10, 16, 16, 10, 14, 20, 20,  2,  2,  2,  2,  2;
        1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1;
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        5,  0, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  5, 12, 12;
        2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  2,  2;
        3,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  1,  1;
    ];

    pdegplot(g,'EdgeLabels','on')
    axis equal

    [p,e,t] = initmesh(g,'hmax',hmax);

    coordinates = p';
    elements    = t(1:3,:)';
    edges       = e(1:2,:)';
    segments    = e(5,:);   % Segmentnummer jeder Randkante

    % Beispiel:
    boundary1Segments = [1 8 14 15 16 17];
    boundary2Segments = [2 3 4 5 6 7 9 10 11 12 13];

    boundary1 = edges(ismember(segments,boundary1Segments),:);
    boundary2 = edges(ismember(segments,boundary2Segments),:);

    boundary = {boundary1,boundary2};
end
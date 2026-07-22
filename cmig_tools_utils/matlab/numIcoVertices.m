function numVertices = numIcoVertices(icoNum)
% Trivial function that gives us how many vertices to pick based on
% icosahedron number
switch icoNum
    case 0
        numVertices = 12;

    case 1
        numVertices = 42;

    case 2
        numVertices = 162;

    case 3
        numVertices = 642;

    case 4
        numVertices = 2562;

    case 5
        numVertices = 10242;

    case 6
        numVertices = 40962;

    case 7
        numVertices = 163842;
end
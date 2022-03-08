function test
    A = [5,4,3,0,-3,-4,-5,-4,-3,0,3,4; 0,3,4,5,4,3,0,-3,-4,-5,-4,-3];
    B = [0,1,2,3,4,5;0,0,0,0,0,0];
    
    altitude(A)
endfunction

function z = altitude (X)
## input format: way-point = [x1, x2, ...; y1, y2, ...]
    p1 = 0.1501;
    p2 = 0.2311;
    p3 = 0.7068;
    p4 = 0.4860;
    p5 = 0.6913;
    p6 = 0.2621;
    p7 = 0.4565;
    x = X(1,:);
    y = X(2,:);
    z = sin(y+p1) + p2.*sin(x) + p3.*cos(p4.*sqrt(x.^2+y.^2)) + p5.*cos(y) + p6.*sin(p6.*sqrt(x.^2+y.^2)) + p7.*cos(y);
endfunction
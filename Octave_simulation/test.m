function test
    % threat_source = [6.2280, 17.781, 15.681, 6.5280, 22.581, 15.057, 21.036; 8.5230, 4.6080, 17.208, 13.629, 21.108, 11.835, 15.846; 2.2826, 1.9663, 2.8540, 2.0762, 1.9393, 2.4483, 2.4404];
    % threat_source = [6.77, 13.57, 10.13, 16.81, 14.04, 23.45; 23.07, 22.65, 17.48, 15.55, 11.22, 8.95; 1.93, 2.47, 1.49, 1.88, 2.52, 2.25];
    threat_source = [-3.1; 3.2; 0.4];

    x = linspace(-4,0,100);
    y = linspace(0,4,100);
    [xx, yy] = meshgrid(x, y);
    z = altitude(xx, yy);
    % for i=1:100
    %     for j=1:100
    %         z(i,j) = altitude([x(i);y(j)]);
    %     endfor
    % endfor 
    figure
    contour(x,y,z);
    plotThreat_source (threat_source);
    hold on;

    
endfunction

function z = altitude (x, y)
    a = 0.4;
    b = 2.6;
    d = 1.3;
    g = 2.4;
%     x = way_point(1,:);
%     y = way_point(2,:);
    z = a.*cos(b.*sqrt(x.^2+y.^2))+d.*sin(g.*sqrt(x.^2+y.^2));
endfunction

function plotThreat_source (threat_source)
    hold on
    for i = 1:columns(threat_source)
        plotCircle(threat_source(:,i));
    endfor
    hold off
endfunction

function plotCircle (threat_source)
    t = linspace(0,2*pi,100)'; 
    circsx = threat_source(3).*cos(t) + threat_source(1); 
    circsy = threat_source(3).*sin(t) + threat_source(2); 
    plot(circsx,circsy, "linewidth", 2); 
endfunction

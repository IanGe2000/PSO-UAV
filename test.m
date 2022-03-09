function test(n, N, G_n)
    % n = 10;         # number of particles in a subgroup
    % N = 10;         # number of way-points
    % G_n = 5;        # number of subgroups
    startpoint = [0; 0];
    endpoint = [30; 30];
    position_range = [min([startpoint, endpoint](2,:)), max([startpoint, endpoint](2,:))];
    xIntervals = linspace(startpoint(1),endpoint(1),N);
    # swarm initialization
    swarm = swarmInit(startpoint, endpoint, position_range, n, N, G_n)
    for i = 1:size(swarm,3)
        figure
        plot(xIntervals, swarm(:,:,i)')
    endfor
endfunction

function swarm = swarmInit (startpoint, endpoint, position_range, n, N, G_n)
## input format: swarm size = n by N by G_n
    swarm = rand(n, N, G_n)*10 - 5 + linspace(startpoint(2),endpoint(2),N);
    swarm(:,1,:) = startpoint(2);
    swarm(:,end,:) = endpoint(2);
    swarm = swarm - (swarm < position_range(1)) .* swarm - (swarm > position_range(2)) .* swarm + (swarm > position_range(2))*position_range(2);
endfunction
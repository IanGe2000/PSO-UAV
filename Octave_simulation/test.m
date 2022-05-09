function test
    xIntervals = [0,1,2,3,4];
    swarm(:,:,1) = [0,0,0,2,3;0,0,0,2,3];
    swarm(:,:,2) = [0,0,0,2,3;0,0,0,2,3];
    swarm2Trajectory(swarm,xIntervals)
endfunction

function trajectory = swarm2Trajectory (swarm, xIntervals)
## convert swarm into trajectories (1 trajectory per page, 1 subgroup of trajectories per book)
    trajectory = zeros(3,columns(swarm),rows(swarm),size(swarm,3));
    for i = 1:size(swarm,3)
        for j = 1:rows(swarm)
            trajectory(:,:,j,i) = particle2Trajectory(swarm(j,:,i), xIntervals);
        endfor
    endfor
endfunction

function trajectory = particle2Trajectory (particle, xIntervals)
## convert particle into trajectory
    trajectory = [xIntervals; particle; altitude([xIntervals; particle])];
endfunction

function z = altitude (way_point)
    p1 = 0.1501;
    p2 = 0.2311;
    p3 = 0.7068;
    p4 = 0.4860;
    p5 = 0.6913;
    p6 = 0.2621;
    p7 = 0.4565;
    x = way_point(1,:);
    y = way_point(2,:);
    z = sin(y+p1) + p2.*sin(x) + p3.*cos(p4.*sqrt(x.^2+y.^2)) + p5.*cos(y) + p6.*sin(p6.*sqrt(x.^2+y.^2)) + p7.*cos(y);
endfunction
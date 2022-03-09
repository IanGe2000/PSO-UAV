## Copyright (C) 2022 iange
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} main
##
## @seealso{}
## @end deftypefn

## Author: iange <iange@IANSOMEN>
## Created: 2022-03-04

% ## termimologies:
% ## row, column, page, book
% ## 
% ## trajectory = [x1, x2, ..., xN; y1, y2, ..., yN; z1, z2, ..., zN]
% ## course = [x1, x2, ..., xN; y1, y2, ..., yN]
% ## xIntervals = [x1, x2, ..., xN]
% ## particle = [y1, y2, ..., yN]
% ## 
% ## subgroup = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yn1, yn2, ..., ynN]
% ##        with each row being a particle
% ## swarm(:,:,i) = subgroup(i)
% ##        with each page being a subgroup
% ## solution = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yG_n1, yG_n2, ..., yG_nN]
% ##        with each row being the best particle from its subgroup
% ## 
% ## N = number of discrete points along the x axis
% ## G_n = number of subgroups
% ## n = number of particles within each subgroup
% ## 
% ## slice = [x1, x2; y1, y2]
% ## segment = [x1, x2, x3; y1, y2, y3; z1, z2, z3]
% ## 
% ## threat_sourse = [T1, T2, ...; r1, r2, ...] = [x1, x2, ...; y1, y2, ...; r1, r2, ...]

function retval = main
    ## Step 1
    n = 10;         # number of particles in a subgroup
    N = 10;         # number of way-points
    G_n = 5;        # number of subgroups
    startpoint = [0; 0];
    endpoint = [30; 30];
    position_range = [min([startpoint, endpoint](2,:)), max([startpoint, endpoint](2,:))];
    xIntervals = linspace(startpoint(1),endpoint(1),N);
    # swarm initialization
    swarm = swarmInit(startpoint, endpoint, position_range, n, N, G_n)
    % for i = 1:size(swarm,3)
    %     figure
    %     plot(xIntervals, swarm(:,:,i)')
    % endfor
    # constraints
    theta_Tmax = 90;
    theta_Cmax = 45;
    threat_source = [6.2280, 17.781, 15.681, 6.5280, 22.581, 15.057, 21.036; 8.5230, 4.6080, 17.208, 13.629, 21.108, 11.835, 15.846; 2.2826, 1.9663, 2.8540, 2.0762, 1.9393, 2.4483, 2.4404];
    threat_source = [30; 0; 1];
    # PSO parameters
    maxgeneration = 1;
    omega = linspace(0.7, 0.4, maxgeneration);
    phi_p = 0.2;    # cognitive coefficient
    phi_g = 0.2;    # social coefficient
    velocity_range = [-6.0, 6.0];
    # objective function weights
    omega_d = 0.3;
    omega_c = 0.5;
    CT = 4.5;
    N_W = N;

    ## Step 2
    generation = 1;
    # Initialize the particle's best position to be the current swarm and the corresponding objective value to be inf
    P_pos = swarm;
    P_obj = inf(rows(swarm),1,size(swarm,3));
    # Initialize the group best particle to be the 1st particle in its subgroup    
    G_pos = swarm(1,:,:);
    G_obj = inf(1,1,size(swarm,3));
    # The solution is apparently the best particle from each group
    solution = shiftdim(G_pos,1)';

    while generation <= maxgeneration:
        # calculate the objective of this iteration of swarm
        P_objective = F(swarm, xIntervals, solution, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
        # refresh P_pos and P_obj
        P_refresh = P_objective < P_obj;
        for i = 1:size(swarm,3)
            for j = 1:rows(swarm)
                if P_refresh(j,1,i)
                    # particle's best is refreshed
                    P_pos(j,:,i) = swarm(j,:,i);
                    P_obj(j,1,i) = P_objective(j,1,i);
                endif
            endfor
        endfor

        ## Step 3
        # refresh G_pos and G_obj
        [G_objective, group_best_index] = min(P_obj);
        G_refresh = G_objective < G_obj;
        for i = 1:size(swarm,3)
            if G_refresh(1,1,i)
                # group best is refreshed
                G_pos(1,:,i) = swarm(group_best_index(1,1,i),:,i);
                G_obj(1,1,i) = G_objective(1,1,i);
            endif
        endfor
        # regenerates the solution from G_pos
        solution = shiftdim(G_pos,1)'
        # check if there is any duplication of particles
        [~, rowindex] = unique(solution, "rows", "stable");
        missingindex = resetIndex(rowindex);
        if ~isempty(missingindex)
            for i = missingindex
                swarm(:,:,i) = swarmInit(startpoint, endpoint, position_range, n, N, 1);
            endfor
            objective = F(swarm, xIntervals, solution, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
            for i = missingindex
                solution(i,:) = swarm(1,:,i);
                P_pos(:,:,i) = swarm(:,:,i);
                P_obj(:,:,i) = F(swarm, xIntervals, solution, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
            endfor
        endif
    endwhile

endfunction

function missingindex = resetIndex (rowindex)
    missingindex = [];
    j = 1;
    for i = 1:rowindex(end)
        if i != rowindex(j)
            missingindex = [missingindex; i];
        else
            j++;
        endif
    endfor
endfunction

function swarm = swarmInit (startpoint, endpoint, position_range, n, N, G_n)
## input format: swarm size = n by N by G_n
    swarm = rand(n, N, G_n)*10 - 5 + linspace(startpoint(2),endpoint(2),N);
    swarm(:,1,:) = startpoint(2);
    swarm(:,end,:) = endpoint(2);
    swarm = swarm - (swarm < position_range(1)) .* swarm - (swarm > position_range(2)) .* swarm + (swarm > position_range(2))*position_range(2);
endfunction

function d = distance (varargin)
## Calculate a bunch of different distance using function overloading
    [data, opt] = parseparams(varargin);
    if columns(data) == 1
    ## input format: X = [X1, X2, ...]
        X = data{:};
        if strcmpi([opt{:}],'recursive') || isempty(opt)
        ## output format: d = [dis(X1,X2), dis(X2,X3), ...]
            d = distance(X,shift(X,-1,2))(1:end-1);
        else
            error("Unexpected option.")
        endif
    elseif columns(data) == 2    
    ## input format: X = [X1, X2, ...]; Y = [Y1, Y2, ...]
        [X, Y] = data{:};
        if strcmpi([opt{:}],'parallel') || isempty(opt)
        ## output format: d = [dis(X1,Y1), dis(X2,Y2), ...]
            if ~size_equal(X, Y)
                error("Dimension does not agree.")
            else
                d = sqrt(sum((X-Y).^2));
            endif
        elseif strcmpi([opt{:}],'mesh')
        ## output format: d = [dis(X1,Y1), dis(X1,Y2), ...; dis(X2,Y1), dis(X2,Y2), ...; ...]
            if rows(X) != rows(Y)
                error("Dimension does not agree.")
            else
                d = zeros(columns(X), columns(Y));
                for i = 1:columns(X)
                    for j = 1:columns(Y)
                        d(i,j) = sqrt(sum((X(:,i)-Y(:,j)).^2));
                    endfor
                endfor
            endif
        elseif strcmpi([opt{:}],'recursive')
        ## input format: course = X; threat_source = Y
        ## output format: d = [dis(Y1,X1X2), dis(Y1,X2X3), ...; dis(Y2,X1X2), dis(Y2,X2X3), ...; ...]
            if rows(X) != 2 || rows(Y) != 2
                error("Only support point to line distance under recursive mode.")
            else
                d = zeros(columns(Y),columns(X)-1);
                for i = 1:columns(Y)
                    for j = 1:columns(X)-1
                        d(i,j) = abs(det([X(1,j)-X(1,j+1), X(2,j)-X(2,j+1); Y(1,i)-X(1,j+1), Y(2,i)-X(2,j+1)]))/distance(X(:,j),X(:,j+1));
                    endfor
                endfor
            endif
        else
            error("Unexpected option.")
        endif
    elseif columns(data) == 3
    ## input format: X = [X1, X2, ...]; A = [A1, A2, ...]; B = [B1, B2, ...]
        [X, A, B] = data{:};
        if strcmpi([opt{:}],'parallel') || isempty(opt)
        ## output format: d = [dis(X1, A1B1), dis(X2, A2B2), ...]
            if ~size_equal(X, A, B)
                error("Dimension does not agree.")
            elseif rows(X) != 2
                error("Only support point to line distance under 3 inputs.")
            else
                d = zeros(1,columns(X));
                for i = 1:columns(X)
                    d(i) = abs(det([A(1,i)-B(1,i), A(2,i)-B(2,i); X(1,i)-B(1,i), X(2,i)-B(2,i)]))/distance(A(:,i),B(:,i));
                endfor
            endif
        elseif strcmpi([opt{:}],'mesh')
        ## output format: d = [dis(X1, A1B1), dis(X1, A2B2), ...; dis(X2, A1B1), dis(X2, A2B2), ...; ...]
            if ~size_equal(A, B)
                error("Dimension does not agree.")
            elseif rows(X) != 2
                error("Only support point to line distance under 3 inputs.")
            else
                d = zeros(columns(X),columns(A));
                for i = 1:columns(X)
                    for j = 1:columns(A)
                        d(i,j) = abs(det([A(1,j)-B(1,j), A(2,j)-B(2,j); X(1,i)-B(1,j), X(2,i)-B(2,j)]))/distance(A(:,j),B(:,j));
                    endfor
                endfor
            endif
        else
            error("Unexpected option.")
        endif
    endif
endfunction

% function output = cos (X, varargin)
% ## overloads the built-in cos function to accept the length of the triangle and calculate the cosin of that angle
%     if columns(varargin) == 0
%         output = cos(X);
%     elseif columns(varargin) == 2
%         pending = (X^2 + varargin{1}^2 - varargin{2}^2)/(2*X*varargin{1});
%         if pending >= -1 && pending <= 1
%             output = pending;
%         else
%             error("Input sides do not form a triangle.")
%         endif
%     endif
% endfunction

function trajectory = particle2Trajectory (particle, xIntervals)
## convert particle into trajectory
    trajectory = [xIntervals; particle; altitude([xIntervals; particle])];
endfunction

function trajectory = course2Trajectory (course)
## convert course into trajectory
    trajectory = [course; altitude(course)];  
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

function pathDistance = F_d (varargin)
## input format: 1 input: trajectory = [x1, x2, ...; y1, y2, ...; z1, z2, ...]
##               2 inputs: particle = [y1, y2, ..., yN]; xIntervals = [x1, x2, ..., xN]
    if columns(varargin) == 1
        trajectory = varargin{1};
    elseif columns(varargin) == 2
        [particle, xIntervals] = varargin{:};
        trajectory = particle2Trajectory (particle, xIntervals);
    endif
    den = distance(trajectory(:,1), trajectory(:,end));
    d = distance(trajectory);
    num = sum(d);
    pathDistance = num/den;
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

function bool = threatConflict (course, threat_source)
## output format: bool = true if the course conflicts with any threat_source
    % disp(size(threat_source))
    d1 = distance(course, threat_source(1:2,:), 'recursive');   ## distance from each threat_source to each slice       ## (T by N-1)
    d2 = distance(course, threat_source(1:2,:), 'mesh');        ## distance from each threat_source to each way-point   ## (N by T)
    d3 = distance(course);                                      ## length of each slice                                 ## (1 by N-1)

    booltable1 = d2<threat_source(3,:);         ## if each way-point is within the range of each threat_source
    booltable1 = (booltable1 | shift(booltable1,-1,1))(1:end-1,:)';     ## if either way-point of each slice is within the range of each threat_source  ## (T by N-1)
    booltable2 = d1>threat_source(3,:)';        ## if the distance from each threat_source to each slice is greater than the source's radius            ## (T by N-1)
    booltable3 = d3.^2>abs(diff(d2.^2))';       ## if c^2>|a^2-b^2|                                                                                     ## (T by N-1)

    booltable = booltable1 | ~booltable1 & ~booltable2 & booltable3;## if each slice conflicts with each threat_source   ## (T by N-1)
    bool = sum(sum(booltable,"native"),"native");               ## if any slice conflicts with any threat_source         ## (1 by 1)

    % disp("threatConflict")
    % d1
    % d2
    % d3    
    % booltable1
    % booltable2
    % booltable3
    % booltable
    % bool
endfunction

function theta = turningAngle (segment)
## output format: theta is the angle formed by the segment in degree
    dif = diff(segment,1,2);
    theta = acosd((dif(1)*dif(4)+dif(2)*dif(5))/(sqrt((dif(1))^2+(dif(2))^2)*sqrt((dif(4))^2+(dif(5))^2)));
endfunction

function bool = maxTurningAngle (trajectory, theta_Tmax)
## input format: theta_Tmax is in degree
## output format: bool = [P1P2P3, P2P3P4, ...](1 by N-2)
    bool = zeros(1,columns(trajectory)-2);
    for i = 1:columns(trajectory)-2
        bool(i) = turningAngle(trajectory(:,i:i+2)) < theta_Tmax;
    endfor

    % disp("maxTurningAngle")
    % bool
endfunction

function cooperativeFactor = F_c (particle, groupindex, solution, CT)
## output format: cooperativeFactor is a nonnegative integer
    D_M = zeros(rows(solution)-1, columns(solution));
    k = 1;
    for i = 1:rows(solution)
        if i == groupindex
            continue
        else
            D_M(k++,:) = abs(particle - solution(i,:));
        endif
    endfor
    cooperativeFactor = sum(min(D_M)<CT);
    % disp("cooperativeFactor")
    % D_M
    % min(D_M)
    % min(D_M)<CT
    % sum(min(D_M)<CT)
endfunction

function altitude = F_a (trajectory, N_W)
## input format: N_W is generally larger than N
## this function is currently unclear about the second parameter: N_W, as how is it larger than N
## and by how much, how to determine the way-points with different intervals along the x Axis.
## As a solution for the time, N_W is treated equal to N.
    altitude = sum(trajectory(3,:))/columns(trajectory);
endfunction

function theta = climbingDivingAngle (slice)
## output format: theta is the angle between the slice and the horison in degree
    theta = atand(abs(diff(altitude(slice)))/sqrt(sum(diff(slice,2).^2)));
endfunction

function bool = maxClimbingDivingAngle (course, theta_Cmax)
## input format: theta_Cmax is in degree
## output format: bool = [P1P2, P2P3, ...](1 by N-1) 
    bool = zeros(1,columns(course)-1);
    for i = 1:columns(course)-1
        bool(i) = climbingDivingAngle(course(:,i:i+1)) < theta_Cmax;
    endfor

    % disp("maxClimbingDivingAngle")
    % bool
endfunction

function objective = F (swarm, xIntervals, solution, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W)
## output format: objective is a G_n by 1 matrix
    trajectory = swarm2Trajectory (swarm, xIntervals);
    objective = zeros(rows(swarm),1, size(swarm,3));
    for i = 1:size(swarm,3)
        for j = 1:rows(swarm)
            if threatConflict([xIntervals; swarm(j,:,i)], threat_source) || ~prod(maxTurningAngle(trajectory(:,:,j,i), theta_Tmax), "native") || ~prod(maxClimbingDivingAngle(trajectory(1:2,:,j,i), theta_Cmax), "native")
                objective(j,1,i) = inf;
            else
                objective(j,1,i) = omega_d * F_d(trajectory(:,:,j,i)) + omega_c * F_c(swarm(j,:,i), i, solution, CT) + (1-omega_c-omega_d) * F_a(trajectory(:,:,j,i), N_W);
            endif
        endfor
    endfor
endfunction

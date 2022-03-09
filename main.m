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

## termimologies:
## trajectory = [x1, x2, ..., xN; y1, y2, ..., yN; z1, z2, ..., zN]
## course = [x1, x2, ..., xN; y1, y2, ..., yN]
## xIntervals = [x1, x2, ..., xN]
## particle = [y1, y2, ..., yN]
## 
## subgroup = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yn1, yn2, ..., ynN]
##        with each row being a particle
## swarm(:,:,i) = subgroup(i)
##        with each page being a subgroup
## solution = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yG_n1, yG_n2, ..., yG_nN]
##        with each row being the best particle from its subgroup
## 
## N = number of discrete points along the x axis
## G_n = number of subgroups
## n = number of particles within each subgroup
## 
## slice = [x1, x2; y1, y2]
## segment = [x1, x2, x3; y1, y2, y3; z1, z2, z3]
## 
## threat_sourse = [T1, T2, ...; r1, r2, ...] = [x1, x2, ...; y1, y2, ...; r1, r2, ...]

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
    swarm = swarmInit(startpoint, endpoint, position_range, n, N, G_n);
    for i = 1:size(swarm,3)
        figure
        plot(xIntervals, swarm(:,:,i)')
    endfor
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
    solution = shiftdim(swarm(1,:,:),1)';
    P_pos = swarm;      # particle's best position
    objective = F(swarm, xIntervals, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
    P_obj = objective;  # particle's best objective
    [G_obj, swarm_best_index] = min(objective); # swarm's best objective
    G_pos = swarm(swarm_best_index,:);          # swarm's best position
    while generation <= maxgeneration
        objective = F(swarm, xIntervals, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
        [objective_min, swarm_best_index] = min(objective);

        generation++;
    endwhile

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
## convert swarm into trajectories (each page is a trajectory)
    trajectory = zeros(3,columns(swarm),rows(swarm));
    for i = 1:rows(swarm)
        trajectory(:,:,i) = particle2Trajectory(swarm(i,:), xIntervals);
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

function bool = threatConflict (course, threat_source)
## input format: course = [P1, P2, ...] = [x1, x2, ...; y1, y2, ...](2 by N), threat_source = [T1, T2, ...; r1, r2, ...] = [x1, x2, ...; y1, y2, ...; r1, r2, ...](3 by T)
## output format: bool (1 by 1) (if the course conflicts with any threat_source)
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
## input format: segment = [x(i-1), x(i), x(i+1); y(i-1), y(i), y(i+1); z(i-1), z(i), z(i+1)]
    dif = diff(segment,1,2);
    theta = acosd((dif(1)*dif(4)+dif(2)*dif(5))/(sqrt((dif(1))^2+(dif(2))^2)*sqrt((dif(4))^2+(dif(5))^2)));
endfunction

function bool = maxTurningAngle (trajectory, theta_Tmax)
## input format: trajectory = [x1, x2, ...; y1, y2, ...; z1, z2, ...]; theta_Tmax is in degree
## output format: bool = [P1P2P3, P2P3P4, ...](1 by N-2)
    bool = zeros(1,columns(trajectory)-2);
    for i = 1:columns(trajectory)-2
        bool(i) = turningAngle(trajectory(:,i:i+2)) < theta_Tmax;
    endfor

    % disp("maxTurningAngle")
    % bool
endfunction

function cooperativeFactor = F_c (swarm, CT, varargin)
## input format: swarm = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yG_n1, yG_n2, ..., yG_nN]
##               CT = cooperative threshold 
##               (particleindex) specifies the particle
## output format: cooperativeFactor is a (G_n by 1) 0/1 matrix when particle is not specified, 0/1 number when it is
    if columns(varargin) == 0
        D_M = zeros(rows(swarm)-1, columns(swarm), rows(swarm));
        k = 1;
        for i = 1:rows(swarm) #page index of D_M / for each particle
            for j = 1:rows(swarm)  #row index of D_M / compared to every other particle
                if j == i
                    continue
                else
                    D_M(k++,:,i) = abs(swarm(i,:) - swarm(j,:));
                endif
            endfor
        endfor
        cooperativeFactor = reshape(sum(min(D_M)>=CT), rows(swarm), 1);
    elseif columns(varargin) == 1
        particleindex = varargin{1};
        D_M = zeros(rows(swarm)-1, columns(swarm));
        k = 1;
        for j = 1:rows(swarm)  #row index of D_M / compared to every other particle
            if j == particleindex
                continue
            else
                D_M(k++,:) = abs(swarm(particleindex,:) - swarm(j,:));
            endif
        endfor
        cooperativeFactor = sum(min(D_M)>=CT);
    endif
    % disp("cooperativeFactor")
    % D_M
    % min(D_M)
    % min(D_M)>=CT
    % sum(min(D_M)>=CT)
endfunction

function altitude = F_a (trajectory, N_W)
## input format: trajectory = [x1, x2, ...; y1, y2, ...; z1, z2, ...]; N_W is generally larger than N
## this function is currently unclear about the second parameter: N_W, as how is it larger than N
## and by how much, how to determine the way-points with different intervals along the x Axis.
## As a solution for the time, N_W is treated equal to N.
    altitude = sum(trajectory(3,:))/columns(trajectory);
endfunction

function theta = climbingDivingAngle (slice)
## input format: slice = [x1, x2;y1, y2]
    theta = atand(abs(diff(altitude(slice)))/sqrt(sum(diff(slice,2).^2)));
endfunction

function bool = maxClimbingDivingAngle (course, theta_Cmax)
## input format: course = [P1, P2, ...] = [x1, x2, ...; y1, y2, ...](2 by N), theta_Cmax is in degree
## output format: bool = [P1P2, P2P3, ...](1 by N-1) 
    bool = zeros(1,columns(course)-1);
    for i = 1:columns(course)-1
        bool(i) = climbingDivingAngle(course(:,i:i+1)) < theta_Cmax;
    endfor

    % disp("maxClimbingDivingAngle")
    % bool
endfunction

function objective = F (swarm, xIntervals, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W)
## input format: swarm = [y11, y12, ..., y1N; y21, y22, ..., y2N; ...; yG_n1, yG_n2, ..., yG_nN]; xIntervals = [x1, x2, ..., xN]
## output format: objective is a G_n by 1 matrix
    trajectory = swarm2Trajectory (swarm, xIntervals);
    objective = zeros(rows(swarm),1);
    for i = 1:rows(swarm)
        % i
        if threatConflict([xIntervals; swarm(i,:)], threat_source) || ~prod(maxTurningAngle(trajectory(:,:,i), theta_Tmax), "native") || ~prod(maxClimbingDivingAngle(trajectory(1:2,:,i), theta_Cmax), "native")
            objective(i) = inf;
        else
            objective(i) = omega_d * F_d(trajectory(:,:,i)) + omega_c * F_c(swarm, CT, i) + (1-omega_c-omega_d) * F_a(trajectory(:,:,i), N_W);
        endif
        % objective(i)
    endfor
endfunction

% Georg Schildbach, 15/Dec/2022 --- RRT+CCD Algorithm
% CCD + RRT path planning method
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2023, Georg Schildbach.
% --------------------------------------------------------------------------------------------------
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
% associated documentation files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or
% substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
% NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% --------------------------------------------------------------------------------------------------
% OUTPUTS:
% - car: structure with car information
% - nObs: number of obstacles (1)
% - sObs: number of vertices of each obstacle (1 x nObs)
% - Obst: array with obstacle vertices in counter-clockwise order (nObs x 2*max(sObs))
% - Nodes: array with path information
%          [x , y , psi , distance to start , parent , distance to target] (nNod x 6)
%          Nodes(1,:): target node
%          Nodes(2,:): initial node
% - nNod: number of nodes, including target and initial node (1)
% - aNod: most advanced node, based on estimated distance to target (1)
% --------------------------------------------------------------------------------------------------

clc
close all
clear all

mex collcheck.c

% 1) Inputs  ---------------------------------------------------------------------------------------

% 1.1) Vehicle

car.w = 2.0;         % width of the car [m]
car.lf = 1.6;        % distance from reference point to front axle [m]
car.lr = 0.9;        % distance from reference point to rear axle [m]
car.lff = 2.1;       % distance from reference point to front of the car [m]
car.lrr = 1.4;       % distance from reference point to rear of the car [m]
car.delta_max = 30;  % maximum front wheel steering angle [deg]
car.R_min = (car.lf+car.lr)/tand(car.delta_max)*cos(atan(car.lr/(car.lf+car.lr)*tand(car.delta_max))); 
                     % minimum turning radius [m]

% 1.2) Map

switch 2
    case 1 % parking scenario
        x_min = 0;    % map x_min [m]
        x_max = 40;   % map x_max [m]
        y_min = 0;    % map y_min [m]
        y_max = 20;   % map y_max [m]
        nObs = 6;     % number of obstacles
        sObs = [4 , 4 , 4 , 4 , 4 , 4];  % number of vertex points for each obstacle
        Obst = [x_min-1 , y_min-1 , x_max+1 , y_min-1 , x_max+1 , y_min , x_min-1 , y_min; % obstacle 1 (lower wall)
                x_min-1 , y_min-1 , x_min , y_min-1 , x_min , y_max+1 , x_min-1 , y_max+1; % obstacle 2 (left wall)
                x_min-1 , y_max , x_max+1 , y_max , x_max+1 , y_max+1 , x_min-1 , y_max+1; % obstacle 3 (upper wall)
                x_max , y_min-1 , x_max+1 , y_min-1 , x_max+1 , y_max+1 , x_max , y_max+1; % obstacle 4 (right wall)
                9.5 , -1 , 10.5 , -1 , 10.5 , 12 , 9.5 , 12;                               % obstacle 5 (left obstacle)
                29.5 , 21 , 29.5 , 8 , 30.5 , 8 , 30.5 , 21];                              % obstacle 6 (right obstacle)
        Init = [3 , 3 , pi/2];                      % initial position: x [m] , y [m] , psi [rad]
        Targ = [35.5 , 16 ,  pi/2];     % target positions: x [m] , y [m] , psi [rad] 
        DeltaTarg = [0.2 , 0.2 , 0.1];  % target interval: Delta x [m] , Delta y [m] , Delta psi [rad]
     case 2 % Scenario 1 (maze)
        x_min = 0;    % map x_min [m]
        x_max = 60;   % map x_max [m]
        y_min = 0;    % map y_min [m]
        y_max = 50;   % map y_max [m]
        nObs = 11;    % number of obstacles
        sObs = linspace(4,4,11);  % number of vertex points for each obstacle
        Obst = [x_min-1 , y_min-1 , x_max+1 , y_min-1 , x_max+1 , y_min , x_min-1 , y_min; % obstacle  1 (lower wall)
                x_min-1 , y_min-1 , x_min , y_min-1 , x_min , y_max+1 , x_min-1 , y_max+1; % obstacle  2 (left wall)
                x_min-1 , y_max , x_max+1 , y_max , x_max+1 , y_max+1 , x_min-1 , y_max+1; % obstacle  3 (upper wall)
                x_max , y_min-1 , x_max+1 , y_min-1 , x_max+1 , y_max+1 , x_max , y_max+1; % obstacle  4 (right wall)
                10 , 10 , 17 , 10 , 17 , 13 , 10 , 13;                                     % obstacle  5 (lower left)
                40 , 10 , 50 , 10 , 50 , 13 , 40 , 13;                                     % obstacle  6 (lower right)
                40 , 20 , 43 , 20 , 43 , 35 , 40 , 35;                                     % obstacle  7 (upper right)
                27 ,  0 , 30 ,  0 , 30 , 25 , 27 , 25;                                     % obstacle  8 (lower wall 1)
                27 , 22 , 27 , 25 , 17 , 25 , 17 , 22;                                     % obstacle  9 (lower wall 2)
                17 , 25 , 20 , 25 , 20 , 37 , 17 , 37;                                     % obstacle 10 (lower wall 3)
                30 , 32 , 33 , 32 , 33 , 50 , 30 , 50];                              % obstacle 11 (upper wall)
        Init = [ 5 ,  5 ,  0];      % initial position: x [m] , y [m] , psi [rad]
        Targ = [52 , 32 ,  0];      % target positions: x [m] , y [m] , psi [rad] 
        DeltaTarg = [0 , 0 , 0.05]; % target interval: Delta x [m] , Delta y [m] , Delta psi [rad]
    case 3
end
    
% 1.3) Options for RRT Algorithm

maxit = 1000;    % maximum number of RRT iterations
maxnod = 500;    % maximum number of constructed nodes

% 2) Modified RRT ----------------------------------------------------------------------------------

rng(1);           % random number seed
Nodes = [Targ , Inf , -1 , 0 ; ...                             % taget node
         Init , 0 , -1 , norm(Targ(1,1:2)-Init(1,1:2),2); ...  % initial node
         zeros(maxnod-1,6) ];  % [x , y , psi , distance from start , parent , distance to goal]
nNod = 2;         % number of nodes
aNod = 2;         % most advanced node (based on distance to target)

% Auxiliary variables

v = zeros(2,1);
q = zeros(2,1);
M = zeros(2,1);
r = zeros(2,1);
phi = zeros(2,1);

for k = 1:maxit
    
    % 2.1) Generate new random node

    r(1) = 100*rand(1);  % random strategy switch
    targetnode = false;

    if Nodes(1,5) == -1  % no connection to target point has been found
        if or( r(1) < 5, k == maxit )  % try target point
            P = Targ(1,1:2)';
            targetnode = true;
        elseif r(1) < 20               % progress from best known node towards target
            r(2) = rand(1);
            P = Nodes(aNod,1:2)' + r(2) * (Targ(1,1:2)' - Nodes(aNod,1:2)');
        elseif r(1) < 40               % progress around most advanced node
            P = [Nodes(aNod,1) + 0.2 * (x_min + (x_max-x_min)*rand(1)) ; ...
                 Nodes(aNod,2) + 0.2 * (y_min + (y_max-y_min)*rand(1))];
        elseif r(1) < 60               % explore around target point
            P = [Targ(1,1) + 0.2 * (x_min + (x_max-x_min)*rand(1)) ; ...
                 Targ(1,2) + 0.2 * (y_min + (y_max-y_min)*rand(1))];
        else                           % explore the whole map
            P = [x_min + (x_max-x_min)*rand(1) ; y_min + (y_max-y_min)*rand(1)];
        end
    else           % a connection to target point has already been found
        if or( r(1) < 10, k == maxit ) % try target point
            P = Targ(1,1:2)';
            targetnode = true;
        elseif r(1) < 50               % explore around target point
            P = [Targ(1,1) + 0.2 * (x_min + (x_max-x_min)*rand(1)) ; ...
                 Targ(1,2) + 0.2 * (y_min + (y_max-y_min)*rand(1))];
        else                           % explore the whole map
            P = [x_min + (x_max-x_min)*rand(1) ; y_min + (y_max-y_min)*rand(1)];
        end
    end

    N = nNod;
    newnode = true;
    
    for i = 2:N  % for every stored node, except the target node
        
        % 2.2) Calculate connecting arc

        % Arc center point M
        
        v = [ cos(Nodes(i,3)) ; sin(Nodes(i,3)) ];
        q = P - Nodes(i,1:2)';  
        if abs( (v(1)*q(2))-(v(2)*q(1)) ) > 10^(-6)  % matrix is non-singular
            M = [ v' ; q' ] \ [ v'*Nodes(i,1:2)' ; 0.5*q'*(P+Nodes(i,1:2)') ];
        else                                         % matrix is singular
            M = (P+Nodes(i,1:2)') / 2 + 10^(+6) * [-q(2);q(1)];
        end

        % Turning direction

        if v' * [ M(2)-Nodes(i,2) ; Nodes(i,1)-M(1) ] >= 0  
            posturn = true;   % vehicle turning counter-clockwise (positive direction)
        else
            posturn = false;  % vehicle turning clockwise (negative direction)
        end
        
        % Arc radius R
        
        R = norm( P - M , 2);
        r(1) = R - car.w/2;  % r_in
        if car.lff >= car.lrr
            r(2) = norm([(R+car.w/2);car.lff],2);  % r_out
        else
            r(2) = norm([(R+car.w/2);car.lrr],2);  % r_out
        end
        
        % Arc limiting angles phi(1), phi(2)
        
        if R >= car.R_min
            if     Nodes(i,1)-M(1) > 0
                phi(1) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1)));
            elseif Nodes(i,1)-M(1) < 0
                phi(1) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1))) + pi;
            else
                if     Nodes(i,2)-M(2) > 0
                    phi(1) = + pi/2;
                elseif Nodes(i,2)-M(2) < 0
                    phi(1) = - pi/2;
                else
                    phi(1) = 0;
                end
            end
            if     P(1)-M(1) > 0
                phi(2) = atan((P(2)-M(2))/(P(1)-M(1)));
            elseif P(1)-M(1) < 0
                phi(2) = atan((P(2)-M(2))/(P(1)-M(1))) + pi;
            else
                if     P(2)-M(2) > 0
                    phi(2) = + pi/2;
                elseif P(2)-M(2) < 0
                    phi(2) = - pi/2;
                else
                    phi(2) = 0;
                end
            end
            if posturn   % vehicle turning counter-clockwise
                phi(2) = phi(2) - floor((phi(2)-phi(1))/2/pi)*2*pi; % phi(2) in [phi(1),phi(1)+2*pi)
            else         % vehicle turning clockwise
                phi(2) = phi(2) - ceil((phi(2)-phi(1))/2/pi)*2*pi;  % phi(2) in (phi(1)-2*pi,phi(1)]
            end
            ang = Nodes(i,3) + (phi(2) - phi(1));
            dist = Nodes(i,4) + R * (phi(2) - phi(1) - floor((phi(2)-phi(1))/2/pi)*2*pi);
        end

        % Check for admissible arc
        
        coll = true;
        if R >= car.R_min
            coll = false;
            if posturn   % vehicle turning counter-clockwise
                phi(1) = phi(1) - asin(car.lrr/(R-car.w/2));
                phi(2) = phi(2) + asin(car.lff/(R-car.w/2));
            else         % vehicle turning clockwise
                q(1) = phi(1);
                phi(1) = phi(2) - asin(car.lff/(R-car.w/2));
                phi(2) = q(1)   + asin(car.lrr/(R-car.w/2));
            end
            for j = 1:nObs
                flag = collcheck(Obst(j,1:2*sObs(j)),M',phi',r');
                if flag
                    coll = true;
                    break
                end
            end
        end
        
        % 2.3) Add node

        % Case 1: Add new node or update node

        if and( ~coll , ~targetnode )
            if or( newnode , dist < Nodes(nNod,4) )
                if newnode
                    nNod = nNod + 1;
                    newnode = false;
                end
                Nodes(nNod,1:6) = [P' , ang , dist , i , norm(Targ(1,1:2)'-P,2)];
                if norm(Targ(1,1:2)'-P,2) < Nodes(aNod,6)
                    aNod = nNod;
                end
            end
        end

        % Case 2: P is target node

        if and( ~coll , targetnode )
            if and( ang>Targ(1,3)-DeltaTarg(1,3) , ang<Targ(1,3)+DeltaTarg(1,3) )
                if dist < Nodes(1,4)  % new optimal path has been found!
                    Nodes(1,4) = dist;
                    Nodes(1,5) = i;
                end
            end
        end

    end
    
    % 2.4) Check for termination
    
    if nNod >= maxnod  % maximum number of nodes reached
        break
    end

end
    

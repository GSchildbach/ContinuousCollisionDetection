% Georg Schildbach, 18/Dec/2022 --- RRT+CCD Algorithm
% Run and plot CCD + RRT path planning method
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
%          [x ; y ; psi ; distance to start ; parent ; distance to target] (6 x (nNod+1))
%          Nodes(:,1): target node
%          Nodes(:,2): initial node
% --------------------------------------------------------------------------------------------------

clc
close all
clear all

mex pathplan.c

% 1) Inputs  ---------------------------------------------------------------------------------------

% 1.1) Vehicle

car = zeros(1,6);
car(1,1) = 2.0;    % width of the car [m]
car(1,2) = 1.6;    % distance from reference point to front axle [m]
car(1,3) = 0.9;    % distance from reference point to rear axle [m]
car(1,4) = 2.1;    % distance from reference point to front of the car [m]
car(1,5) = 1.4;    % distance from reference point to rear of the car [m]
car(1,6) = 30;     % maximum front wheel steering angle [deg]

% 1.2) Map

switch 3
    case 1 % parking scenario
        xlim = zeros(1,4);
        xlim(1,1) = 0;    % map x_min [m]
        xlim(1,2) = 40;   % map x_max [m]
        xlim(1,3) = 0;    % map y_min [m]
        xlim(1,4) = 20;   % map y_max [m]
        nObs = 6;     % number of obstacles
        sObs = [4 , 4 , 4 , 4 , 4 , 4];  % number of vertex points for each obstacle
        Obst = [xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3) , xlim(1,1)-1 , xlim(1,3), ... % obstacle 1 (lower wall)
                xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,1) , xlim(1,3)-1 , xlim(1,1) , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle 2 (left wall)
                xlim(1,1)-1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle 3 (upper wall)
                xlim(1,2) , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,2) , xlim(1,4)+1, ... % obstacle 4 (right wall)
                9.5 , -1 , 10.5 , -1 , 10.5 , 12 , 9.5 , 12, ...                               % obstacle 5 (left obstacle)
                29.5 , 21 , 29.5 , 8 , 30.5 , 8 , 30.5 , 21];                                  % obstacle 6 (right obstacle)
        Init = [3 , 3 , pi/2];                      % initial position: x [m] , y [m] , psi [rad]
        Targ = [35.5 , 16 ,  pi/2];     % target positions: x [m] , y [m] , psi [rad] 
        DeltaTarg = [0.2 , 0.2 , 0.1];  % target interval: Delta x [m] , Delta y [m] , Delta psi [rad]
     case 2 % Scenario 1 (Maze)
        xlim = zeros(1,4);
        xlim(1,1) = 0;    % map x_min [m]
        xlim(1,2) = 60;   % map x_max [m]
        xlim(1,3) = 0;    % map y_min [m]
        xlim(1,4) = 50;   % map y_max [m]
        nObs = 11;    % number of obstacles
        sObs = linspace(4,4,11);  % number of vertex points for each obstacle
        Obst = [xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3) , xlim(1,1)-1 , xlim(1,3), ... % obstacle  1 (lower wall)
                xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,1) , xlim(1,3)-1 , xlim(1,1) , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle  2 (left wall)
                xlim(1,1)-1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle  3 (upper wall)
                xlim(1,2) , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,2) , xlim(1,4)+1, ... % obstacle  4 (right wall)
                10 , 10 , 17 , 10 , 17 , 13 , 10 , 13, ...                                     % obstacle  5 (lower left)
                40 , 10 , 50 , 10 , 50 , 13 , 40 , 13, ...                                     % obstacle  6 (lower right)
                40 , 20 , 43 , 20 , 43 , 35 , 40 , 35, ...                                     % obstacle  7 (upper right)
                27 ,  0 , 30 ,  0 , 30 , 25 , 27 , 25, ...                                     % obstacle  8 (lower wall 1)
                27 , 22 , 27 , 25 , 17 , 25 , 17 , 22, ...                                     % obstacle  9 (lower wall 2)
                17 , 25 , 20 , 25 , 20 , 37 , 17 , 37, ...                                     % obstacle 10 (lower wall 3)
                30 , 32 , 33 , 32 , 33 , 50 , 30 , 50];                              % obstacle 11 (upper wall)
        Init = [ 5 ,  5 ,  0];      % initial position: x [m] , y [m] , psi [rad]
        Targ = [52 , 32 ,  0];      % target positions: x [m] , y [m] , psi [rad] 
        DeltaTarg = [0 , 0 , 0.05]; % target interval: Delta x [m] , Delta y [m] , Delta psi [rad]
    case 3 % Scenario 2 (Obstacles)
        xlim = zeros(1,4);
        xlim(1,1) = 0;    % map x_min [m]
        xlim(1,2) = 60;   % map x_max [m]
        xlim(1,3) = 0;    % map y_min [m]
        xlim(1,4) = 50;   % map y_max [m]
        nObs = 23;    % number of obstacles
        sObs = linspace(4,4,23);  % number of vertex points for each obstacle
        Obst = [xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3) , xlim(1,1)-1 , xlim(1,3), ... % obstacle  1 (lower wall)
                xlim(1,1)-1 , xlim(1,3)-1 , xlim(1,1) , xlim(1,3)-1 , xlim(1,1) , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle  2 (left wall)
                xlim(1,1)-1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4) , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,1)-1 , xlim(1,4)+1, ... % obstacle  3 (upper wall)
                xlim(1,2) , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,3)-1 , xlim(1,2)+1 , xlim(1,4)+1 , xlim(1,2) , xlim(1,4)+1, ... % obstacle  4 (right wall)
                22.5 , 42.5 , 25 , 42.5 , 25 , 45 , 22.5 , 45, ...      % obstacle  1 (from top left...)
                10 , 37.5 , 12.5 , 37.5 , 12.5 , 40 , 10 , 40, ...      % obstacle  2 
                33.5 , 37.5 , 36 , 37.5 , 36 , 40 , 33.5 , 40, ...      % obstacle  3
                26.5 , 35 , 29 , 35 , 29 , 37.5 , 26.5 , 37.5, ...      % obstacle  4
                20 , 32.5 , 22.5 , 32.5 , 22.5 , 35 , 20 , 35, ...      % obstacle  5
                42.5 , 32.5 , 45 , 32.5 , 45 , 35 , 42.5 , 35, ...      % obstacle  6
                31.5 , 31.5 , 34 , 31.5 , 34 , 34 , 31.5 , 34, ...      % obstacle  7
                23.5 , 27.5 , 26 , 27.5 , 26 , 30 , 23.5 , 30, ...      % obstacle  8
                37.5 , 26.5 , 40 , 26.5 , 40 , 29 , 37.5 , 29, ...      % obstacle  9
                46.5 , 26.5 , 49 , 26.5 , 49 , 29 , 46.5 , 29, ...      % obstacle 10
                30 , 25 , 32.5 , 25 , 32.5 , 27.5 , 30 , 27.5, ...      % obstacle 11
                19 , 23 , 21.5 , 23 , 21.5 , 25.5 , 19 , 25.5, ...      % obstacle 12
                25 , 21 , 27.5 , 21 , 27.5 , 23.5 , 25 , 23.5, ...      % obstacle 13
                34 , 20 , 36.5 , 20 , 36.5 , 22.5 , 34 , 22.5, ...      % obstacle 14
                42.5 , 20 , 45 , 20 , 45 , 22.5 , 42.5 , 22.5, ...      % obstacle 15
                27.5 , 15 , 30 , 15 , 30 , 17.5 , 27.5 , 17.5, ...      % obstacle 16
                35 , 15 , 37.5 , 15 , 37.5 , 17.5 , 35 , 17.5, ...      % obstacle 17
                52.5 , 12.5 , 55 , 12.5 , 55 , 15 , 52.5 , 15, ...      % obstacle 18
                40 , 10 , 42.5 , 10 , 42.5 , 12.5 , 40 , 12.5];         % obstacle 19 (...to bottom right)
        Init = [ 5 ,  5 ,  0];      % initial position: x [m] , y [m] , psi [rad]
        Targ = [52 , 38 ,  0];      % target positions: x [m] , y [m] , psi [rad] 
        DeltaTarg = [0 , 0 , 0.05]; % target interval: Delta x [m] , Delta y [m] , Delta psi [rad]
end

% 2) Formats and Settings --------------------------------------------------------------------------

% 2.1) Final path

plotpath.linestyle = '-';
plotpath.linewidth = 2.0;
plotpath.color = [1 0 0]; % RGB value
plotpath.marker = '.';
plotpath.markersize = 6;
plotpath.markeredgecolor = [0 0 0];
plotpath.markerfacecolor = [0 0 0];

% 2.2) Planned paths

plotplan.linestyle = '-';
plotplan.linewidth = 0.5;
plotplan.color = [0.5 0.5 0.5]; % RGB value
plotplan.marker = '.';
plotplan.markersize = 4;
plotplan.markeredgecolor = [0.5 0.5 0.5];
plotplan.markerfacecolor = [0.5 0.5 0.5];

% 2.3) Initial node and target node

plotnode.linestyle = '-';
plotnode.linewidth = 2;
plotnode.color = [0 0 1]; % RGB value
plotnode.marker = '.';
plotnode.markersize = 10;
plotnode.markeredgecolor = [0 0 1];
plotnode.markerfacecolor = [0 0 1];

% 3) Run Path Planner  -----------------------------------------------------------------------------

tic
Nodes = pathplan(car,xlim,sObs,Obst,Init,Targ,DeltaTarg);
T = toc;

display(['Total time for path planning: ' num2str(T,'%2.2f') ' seconds.'])

% 4) Plot Map --------------------------------------------------------------------------------------

% 4.1) Initialize Plot

fig = figure()
hold on
axis equal
xlabel('x [m]')
ylabel('y [m]')
title('RRT Path (gray: planned paths, red: optimal path)')

% 4.2) Find Plot Range

xyrange = [+Inf, -Inf, +Inf, -Inf];
for k = 1:nObs
    2*sum(sObs(1,1:(k-1)));

    if min(Obst(1,2*sum(sObs(1,1:(k-1)))+1:2:2*sum(sObs(1,1:k)))') < xyrange(1)
        xyrange(1) = min(Obst(1,2*sum(sObs(1,1:(k-1)))+1:2:2*sum(sObs(1,1:k)))');
    end
    if max(Obst(1,2*sum(sObs(1,1:(k-1)))+1:2:2*sum(sObs(1,1:k)))') > xyrange(2)
        xyrange(2) = max(Obst(1,2*sum(sObs(1,1:(k-1)))+1:2:2*sum(sObs(1,1:k)))');
    end
    if min(Obst(1,2*sum(sObs(1,1:(k-1)))+2:2:2*sum(sObs(1,1:k)))') < xyrange(3)
        xyrange(3) = min(Obst(1,2*sum(sObs(1,1:(k-1)))+2:2:2*sum(sObs(1,1:k)))');
    end
    if max(Obst(1,2*sum(sObs(1,1:(k-1)))+2:2:2*sum(sObs(1,1:k)))') > xyrange(4)
        xyrange(4) = max(Obst(1,2*sum(sObs(1,1:(k-1)))+2:2:2*sum(sObs(1,1:k)))');
    end
end 
if (xyrange(1,2)-xyrange(1,1))>=(xyrange(1,4)-xyrange(1,3))
    x1 = xyrange(1,1);
    x2 = xyrange(1,2);
    y1 = xyrange(1,3)-(xyrange(1,2)-xyrange(1,1)-xyrange(1,4)+xyrange(1,3))/2;
    y2 = xyrange(1,4)+(xyrange(1,2)-xyrange(1,1)-xyrange(1,4)+xyrange(1,3))/2;
else
    x1 = xyrange(1,1)-(xyrange(1,4)-xyrange(1,3)-xyrange(1,2)+xyrange(1,1))/2;
    x2 = xyrange(1,2)+(xyrange(1,4)-xyrange(1,3)-xyrange(1,2)+xyrange(1,1))/2;
    y1 = xyrange(1,3);
    y2 = xyrange(1,4);
end

if 1  % manual adjustments
    x1 = x1 - 1;
    x2 = x2 + 1;
    y1 = y1 + 4;
    y2 = y2 - 4;
end

% 4.3) Plot Boundaries and Obstacles

fill([x1,xyrange(1,1),xyrange(1,1),x1],[y1,y1,y2,y2],[0 0 0],'EdgeColor',[0 0 0]);
fill([x2,xyrange(1,2),xyrange(1,2),x2],[y1,y1,y2,y2],[0 0 0],'EdgeColor',[0 0 0]);
fill([x1,x2,x2,x1],[y1,y1,xyrange(1,3),xyrange(1,3)],[0 0 0],'EdgeColor',[0 0 0]);
fill([x1,x2,x2,x1],[y2,y2,xyrange(1,4),xyrange(1,4)],[0 0 0],'EdgeColor',[0 0 0]);
for k = 1:nObs
    fill(Obst(1,2*sum(sObs(1,1:(k-1)))+1:2:2*sum(sObs(1,1:k))), ...
         Obst(1,2*sum(sObs(1,1:(k-1)))+2:2:2*sum(sObs(1,1:k))),[0 0 0],'EdgeColor',[0 0 0]);
end
axis([x1,x2,y1,y2]);

% 5) Plot Paths --------------------------------------------------------------------------------

% 5.1) Find nNod

nNod = 2;

while nNod < size(Nodes,2)
    if Nodes(5,nNod+1) ~=0
        nNod = nNod + 1;
    else
        break
    end
end

% 5.2) Plot All Planned Paths

for i = 3:nNod

    % Arc center point M

    v = [ cos(Nodes(3,Nodes(5,i))) ; sin(Nodes(3,Nodes(5,i))) ];
    q = Nodes(1:2,i) - Nodes(1:2,Nodes(5,i));
    if abs( (v(1)*q(2))-(v(2)*q(1)) ) > 10^(-6)  % matrix is non-singular
        M = [ v' ; q' ] \ [ v'*Nodes(1:2,Nodes(5,i)) ; 0.5*q'*(Nodes(1:2,i)+Nodes(1:2,Nodes(5,i))) ];
    else                                         % matrix is singular
        M = (Nodes(i,1:2)'+Nodes(1:2,Nodes(5,i))) / 2 + 10^(+6) * [-q(2);q(1)];
    end

    % Turning direction

    if v' * [ M(2)-Nodes(2,Nodes(5,i)) ; Nodes(1,Nodes(5,i))-M(1) ] >= 0  
        posturn = true;   % vehicle turning counter-clockwise (positive direction)
    else
        posturn = false;  % vehicle turning clockwise (negative direction)
    end

    % Arc radius R

    R = norm( Nodes(1:2,i) - M , 2);

    % Arc angles

    if     Nodes(1,Nodes(5,i))-M(1) > 0
        phi(1) = atan((Nodes(2,Nodes(5,i))-M(2))/(Nodes(1,Nodes(5,i))-M(1)));
    elseif Nodes(1,Nodes(5,i))-M(1) < 0
        phi(1) = atan((Nodes(2,Nodes(5,i))-M(2))/(Nodes(1,Nodes(5,i))-M(1))) + pi;
    else
        if     Nodes(2,Nodes(5,i))-M(2) > 0
            phi(1) = + pi/2;
        elseif Nodes(2,Nodes(5,i))-M(2) < 0
            phi(1) = - pi/2;
        else
            phi(1) = 0;
        end
    end
    if     Nodes(1,i)-M(1) > 0
        phi(2) = atan((Nodes(2,i)-M(2))/(Nodes(1,i)-M(1)));
    elseif Nodes(1,i)-M(1) < 0
        phi(2) = atan((Nodes(2,i)-M(2))/(Nodes(1,i)-M(1))) + pi;
    else
        if     Nodes(2,i)-M(2) > 0
            phi(2) = + pi/2;
        elseif Nodes(2,i)-M(2) < 0
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
    
    % Plot planned paths

    str = plotplan.marker;
    plotplan.marker = 'none';
    cosphi = cos( phi(1) + linspace(0,1,21) * (phi(2)-phi(1)) );
    sinphi = sin( phi(1) + linspace(0,1,21) * (phi(2)-phi(1)) );
    plot([M(1) + R * cosphi] , [M(2) + R * sinphi] , plotplan);
    plotplan.marker = str;
    str = plotplan.linestyle;
    plotplan.linestyle = 'none';
    plot([M(1) + R * cos(phi(1)) , M(1) + R * cos(phi(2))], ...
         [M(2) + R * sin(phi(1)) , M(2) + R * sin(phi(2))],plotplan);
    plotplan.linestyle = str;

end

% 3.2) Plot Optimal Path

% Initial node and target node

plot(Init(1),Init(2),plotnode);
plot(Targ(1),Targ(2),plotnode);
quiver(Init(1),Init(2),(car(4)+car(5))*cos(Init(3)),(car(4)+car(5))*sin(Init(3)),'LineWidth',plotnode.linewidth,'Color',plotnode.color,'MarkerSize',plotnode.markersize)
quiver(Targ(1),Targ(2),(car(4)+car(5))*cos(Targ(3)),(car(4)+car(5))*sin(Targ(3)),'LineWidth',plotnode.linewidth,'Color',plotnode.color,'MarkerSize',plotnode.markersize)
%fill([Targ(1,1),Targ(2,1),Targ(2,1),Targ(1,1)],[Targ(1,2),Targ(1,2),Targ(2,2),Targ(2,2)],plotnode.color,'EdgeColor',plotnode.color);

% Optimal path
  
i = 1;    % start with target node

while Nodes(5,i) > 0

    % Arc center point M

    v = [ cos(Nodes(3,Nodes(5,i))) ; sin(Nodes(3,Nodes(5,i))) ];
    q = Nodes(1:2,i) - Nodes(1:2,Nodes(5,i));
    if abs( (v(1)*q(2))-(v(2)*q(1)) ) > 10^(-6)  % matrix is non-singular
        M = [ v' ; q' ] \ [ v'*Nodes(1:2,Nodes(5,i)) ; 0.5*q'*(Nodes(1:2,i)+Nodes(1:2,Nodes(5,i))) ];
    else                                         % matrix is singular
        M = (Nodes(1:2,i)+Nodes(1:2,Nodes(5,i))) / 2 + 10^(+6) * [-q(2);q(1)];
    end

    % Turning direction

    if v' * [ M(2)-Nodes(2,Nodes(5,i)) ; Nodes(1,Nodes(5,i))-M(1) ] >= 0  
        posturn = true;   % vehicle turning counter-clockwise (positive direction)
    else
        posturn = false;  % vehicle turning clockwise (negative direction)
    end

    % Arc radius R

    R = norm( Nodes(1:2,i) - M , 2);

    % Arc angles

    if     Nodes(1,Nodes(5,i))-M(1) > 0
        phi(1) = atan((Nodes(2,Nodes(5,i))-M(2))/(Nodes(1,Nodes(5,i))-M(1)));
    elseif Nodes(1,Nodes(5,i))-M(1) < 0
        phi(1) = atan((Nodes(2,Nodes(5,i))-M(2))/(Nodes(1,Nodes(5,i))-M(1))) + pi;
    else
        if     Nodes(2,Nodes(5,i))-M(2) > 0
            phi(1) = + pi/2;
        elseif Nodes(2,Nodes(5,i))-M(2) < 0
            phi(1) = - pi/2;
        else
            phi(1) = 0;
        end
    end
    if     Nodes(1,i)-M(1) > 0
        phi(2) = atan((Nodes(2,i)-M(2))/(Nodes(1,i)-M(1)));
    elseif Nodes(1,i)-M(1) < 0
        phi(2) = atan((Nodes(2,i)-M(2))/(Nodes(1,i)-M(1))) + pi;
    else
        if     Nodes(2,i)-M(2) > 0
            phi(2) = + pi/2;
        elseif Nodes(2,i)-M(2) < 0
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
    
    % Plot optimal path

    str = plotpath.marker;
    plotpath.marker = 'none';
    cosphi = cos( phi(1) + linspace(0,1,21) * (phi(2)-phi(1)) );
    sinphi = sin( phi(1) + linspace(0,1,21) * (phi(2)-phi(1)) );
    plot([M(1) + R * cosphi] , [M(2) + R * sinphi] , plotpath);
    plotpath.marker = str;
    str = plotpath.linestyle;
    plotpath.linestyle = 'none';
    plot([M(1) + R * cos(phi(1)) , M(1) + R * cos(phi(2))], ...
         [M(2) + R * sin(phi(1)) , M(2) + R * sin(phi(2))],plotpath);
    plotpath.linestyle = str;
    
    i = Nodes(5,i);

end

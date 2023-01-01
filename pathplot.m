% Georg Schildbach, 17/Dec/2022 --- RRT+CCD Algorithm
% Plots the results of the RRT path planning
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
% INPUTS
% - car: structure with car information
% - nObs: number of obstacles (1)
% - sObs: number of vertices of each obstacle (1 x nObs)
% - Obst: array with obstacle vertices in counter-clockwise order (nObs x 2*max(sObs))
% - Nodes: array with path information
%          [x , y , psi , distance to start , parent , distance to target] ((nNod+1) x 6)
%          Nodes(1,:): target node
%          Nodes(2,:): initial node
% - nNod: number of nodes, including target and initial node (1)
% --------------------------------------------------------------------------------------------------
% OUTPUTS
% - plot of all planned paths and the optimal path
% --------------------------------------------------------------------------------------------------

% 1) Formats and Settings --------------------------------------------------------------------------

% 1.1) Final path

plotpath.linestyle = '-';
plotpath.linewidth = 1.5;
plotpath.color = [1 0 0]; % RGB value
plotpath.marker = '.';
plotpath.markersize = 5;
plotpath.markeredgecolor = [0 0 0];
plotpath.markerfacecolor = [0 0 0];

% 1.2) Planned paths

plotplan.linestyle = '-';
plotplan.linewidth = 0.5;
plotplan.color = [0.5 0.5 0.5]; % RGB value
plotplan.marker = '.';
plotplan.markersize = 5;
plotplan.markeredgecolor = [0.5 0.5 0.5];
plotplan.markerfacecolor = [0.5 0.5 0.5];

% 1.3) Initial node and target node

plotnode.linestyle = '-';
plotnode.linewidth = 2;
plotnode.color = [0 0 1]; % RGB value
plotnode.marker = '.';
plotnode.markersize = 10;
plotnode.markeredgecolor = [0 0 1];
plotnode.markerfacecolor = [0 0 1];

% 2) Plot Map --------------------------------------------------------------------------------------

% 2.1) Initialize Plot

fig = figure()
hold on
axis equal
xlabel('x [m]')
ylabel('y [m]')
title('RRT Path (gray: planned paths, red: optimal path)')

% 2.2) Find Plot Range

xyrange = [+Inf, -Inf, +Inf, -Inf];
for k = 1:nObs
    if min(Obst(k,1:2:2*sObs(1,k))') < xyrange(1)
        xyrange(1) = min(Obst(k,1:2:2*sObs(1,k))');
    end
    if max(Obst(k,1:2:2*sObs(1,k))') > xyrange(2)
        xyrange(2) = max(Obst(k,1:2:2*sObs(1,k))');
    end
    if min(Obst(k,2:2:2*sObs(1,k))') < xyrange(3)
        xyrange(3) = min(Obst(k,2:2:2*sObs(1,k))');
    end
    if max(Obst(k,2:2:2*sObs(1,k))') > xyrange(4)
        xyrange(4) = max(Obst(k,2:2:2*sObs(1,k))');
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

if 0  % manual adjustments
    x1 = x1 - 1;
    x2 = x2 + 1;
    y1 = y1 + 4;
    y2 = y2 - 4;
end

% 2.3) Plot Boundaries and Obstacles

fill([x1,xyrange(1,1),xyrange(1,1),x1],[y1,y1,y2,y2],[0 0 0],'EdgeColor',[0 0 0]);
fill([x2,xyrange(1,2),xyrange(1,2),x2],[y1,y1,y2,y2],[0 0 0],'EdgeColor',[0 0 0]);
fill([x1,x2,x2,x1],[y1,y1,xyrange(1,3),xyrange(1,3)],[0 0 0],'EdgeColor',[0 0 0]);
fill([x1,x2,x2,x1],[y2,y2,xyrange(1,4),xyrange(1,4)],[0 0 0],'EdgeColor',[0 0 0]);
for k = 1:nObs
    fill(Obst(k,1:2:2*sObs(1,k)),Obst(k,2:2:2*sObs(1,k)),[0 0 0],'EdgeColor',[0 0 0]);
end
axis([x1,x2,y1,y2]);

% 3) Plot Paths --------------------------------------------------------------------------------

% 3.1) Plot All Planned Paths

for i = 3:nNod

    % Arc center point M

    v = [ cos(Nodes(Nodes(i,5),3)) ; sin(Nodes(Nodes(i,5),3)) ];
    q = Nodes(i,1:2)' - Nodes(Nodes(i,5),1:2)';
    if abs( (v(1)*q(2))-(v(2)*q(1)) ) > 10^(-6)  % matrix is non-singular
        M = [ v' ; q' ] \ [ v'*Nodes(Nodes(i,5),1:2)' ; 0.5*q'*(Nodes(i,1:2)'+Nodes(Nodes(i,5),1:2)') ];
    else                                         % matrix is singular
        M = (Nodes(i,1:2)'+Nodes(Nodes(i,5),1:2)') / 2 + 10^(+6) * [-q(2);q(1)];
    end

    % Turning direction

    if v' * [ M(2)-Nodes(Nodes(i,5),2) ; Nodes(Nodes(i,5),1)-M(1) ] >= 0  
        posturn = true;   % vehicle turning counter-clockwise (positive direction)
    else
        posturn = false;  % vehicle turning clockwise (negative direction)
    end

    % Arc radius R

    R = norm( Nodes(i,1:2)' - M , 2);

    % Arc angles

    if     Nodes(Nodes(i,5),1)-M(1) > 0
        phi(1) = atan((Nodes(Nodes(i,5),2)-M(2))/(Nodes(Nodes(i,5),1)-M(1)));
    elseif Nodes(Nodes(i,5),1)-M(1) < 0
        phi(1) = atan((Nodes(Nodes(i,5),2)-M(2))/(Nodes(Nodes(i,5),1)-M(1))) + pi;
    else
        if     Nodes(Nodes(i,5),2)-M(2) > 0
            phi(1) = + pi/2;
        elseif Nodes(Nodes(i,5),2)-M(2) < 0
            phi(1) = - pi/2;
        else
            phi(1) = 0;
        end
    end
    if     Nodes(i,1)-M(1) > 0
        phi(2) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1)));
    elseif Nodes(i,1)-M(1) < 0
        phi(2) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1))) + pi;
    else
        if     Nodes(i,2)-M(2) > 0
            phi(2) = + pi/2;
        elseif Nodes(i,2)-M(2) < 0
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
    cosphi = cos( phi(1) + linspace(0,1,101) * (phi(2)-phi(1)) );
    sinphi = sin( phi(1) + linspace(0,1,101) * (phi(2)-phi(1)) );
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

while Nodes(i,5) > 0

    % Arc center point M

    v = [ cos(Nodes(Nodes(i,5),3)) ; sin(Nodes(Nodes(i,5),3)) ];
    q = Nodes(i,1:2)' - Nodes(Nodes(i,5),1:2)';
    if abs( (v(1)*q(2))-(v(2)*q(1)) ) > 10^(-6)  % matrix is non-singular
        M = [ v' ; q' ] \ [ v'*Nodes(Nodes(i,5),1:2)' ; 0.5*q'*(Nodes(i,1:2)'+Nodes(Nodes(i,5),1:2)') ];
    else                                         % matrix is singular
        M = (Nodes(i,1:2)'+Nodes(Nodes(i,5),1:2)') / 2 + 10^(+6) * [-q(2);q(1)];
    end

    % Turning direction

    if v' * [ M(2)-Nodes(Nodes(i,5),2) ; Nodes(Nodes(i,5),1)-M(1) ] >= 0  
        posturn = true;   % vehicle turning counter-clockwise (positive direction)
    else
        posturn = false;  % vehicle turning clockwise (negative direction)
    end

    % Arc radius R

    R = norm( Nodes(i,1:2)' - M , 2);

    % Arc angles

    if     Nodes(Nodes(i,5),1)-M(1) > 0
        phi(1) = atan((Nodes(Nodes(i,5),2)-M(2))/(Nodes(Nodes(i,5),1)-M(1)));
    elseif Nodes(Nodes(i,5),1)-M(1) < 0
        phi(1) = atan((Nodes(Nodes(i,5),2)-M(2))/(Nodes(Nodes(i,5),1)-M(1))) + pi;
    else
        if     Nodes(Nodes(i,5),2)-M(2) > 0
            phi(1) = + pi/2;
        elseif Nodes(Nodes(i,5),2)-M(2) < 0
            phi(1) = - pi/2;
        else
            phi(1) = 0;
        end
    end
    if     Nodes(i,1)-M(1) > 0
        phi(2) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1)));
    elseif Nodes(i,1)-M(1) < 0
        phi(2) = atan((Nodes(i,2)-M(2))/(Nodes(i,1)-M(1))) + pi;
    else
        if     Nodes(i,2)-M(2) > 0
            phi(2) = + pi/2;
        elseif Nodes(i,2)-M(2) < 0
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
    cosphi = cos( phi(1) + linspace(0,1,101) * (phi(2)-phi(1)) );
    sinphi = sin( phi(1) + linspace(0,1,101) * (phi(2)-phi(1)) );
    plot([M(1) + R * cosphi] , [M(2) + R * sinphi] , plotpath);
    plotpath.marker = str;
    str = plotpath.linestyle;
    plotpath.linestyle = 'none';
    plot([M(1) + R * cos(phi(1)) , M(1) + R * cos(phi(2))], ...
         [M(2) + R * sin(phi(1)) , M(2) + R * sin(phi(2))],plotpath);
    plotpath.linestyle = str;
    
    i = Nodes(i,5);

end
    
    

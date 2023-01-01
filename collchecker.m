% Georg Schildbach, 21/Apr/2021 --- RRT+CCD Algorithm
% CCD algorithm
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
% INPUTS:
% Obs: List of vertices of the polytope (nObs x 2)
% M: center point of the arc [x;y] (2 x 1)
% phi: initial and final angle of the arc [phi_1;phi_2] (2 x 1)
% r: inner and outer radius of the arc [r_in;r_out] (2 x 1)
% --------------------------------------------------------------------------------------------------
% OUTPUTS:
% coll: collision indicator (= false if no collision, = true if collision)
% --------------------------------------------------------------------------------------------------

function [coll,term] = collchecker(Obs,M,phi,r)

% 1) Preprocessing  --------------------------------------------------------------------------------

n_Obs = size(Obs,1);
if phi(2) < phi(1)
    phi = [phi(2);phi(1)];
end

% 2) Collision Checking  ---------------------------------------------------------------------------

term = false;
coll = false;          % collision indicator variable
IntsP1 = zeros(2,2);   % list of intersection points with line 1
n_IntsP1 = uint8(0);   % number of intersection points with line 1
r_IntsP1 = zeros(1,2); % radius of intersection points with line 1
IntsP2 = zeros(2,2);   % list of intersection points with line 2
n_IntsP2 = uint8(0);   % number of intersection points with line 2
r_IntsP2 = zeros(1,2); % radius of intersection points with line 2

for k = 1:n_Obs  % obstacle edge number k
    v1 = Obs(k  ,1:2)';      % vertex 1
    if k < n_Obs
        v2 = Obs(k+1,1:2)';  % vertex 2
    else
        v2 = Obs(1,1:2)';    % vertex 2
    end

% 2.1) Intersection points

    ints = false;  % intersection indicator variable
    e = v2 - v1;
    d = M - v1;
    if (d(2)*e(1) - d(1)*e(2)) < 0  % switch vertices s.t. v1 has a smaller angle
        w = v1;
        v1 = v2;
        v2 = w;
        e = v2 - v1;
        d = M - v1;
    end
    
    % Cone line 1
    w = [cos(phi(1)) ; sin(phi(1))];  % ab_check = [e, -w] \ r; 
    if     e(1) ~= 0
        beta = (d(2) - e(2)/e(1)*d(1)) / (e(2)/e(1)*w(1) - w(2));
        alpha = d(1)/e(1) + w(1)/e(1)*beta;
    elseif e(2) ~= 0
        beta = (d(1) - e(1)/e(2)*d(2)) / (e(1)/e(2)*w(2) - w(1));
        alpha = d(2)/e(2) + w(2)/e(2)*beta;
    end
    if and(alpha>=0,alpha<=1)  % intersection is within edge
        n_IntsP1 = n_IntsP1 + 1;
        IntsP1(1:2,n_IntsP1) = v1 + alpha*e - M;     % intersection point relative to M
        if beta >= 0
            r_IntsP1(1,n_IntsP1) = + norm(IntsP1(1:2,n_IntsP1),2); % radius of intersection point
        else
            r_IntsP1(1,n_IntsP1) = - norm(IntsP1(1:2,n_IntsP1),2); % radius of intersection point
        end
        if and(r_IntsP1(1,n_IntsP1)>=r(1),r_IntsP1(1,n_IntsP1)<=r(2))
            ints = true;
            coll = true;
        end
    end
    
    % Cone line 2
    w = [cos(phi(2)) ; sin(phi(2))];  % ab_check = [e, -w] \ d; 
    if     e(1) ~= 0
        beta = (d(2) - e(2)/e(1)*d(1)) / (e(2)/e(1)*w(1) - w(2));
        alpha = d(1)/e(1) + w(1)/e(1)*beta;
    elseif e(2) ~= 0
        beta = (d(1) - e(1)/e(2)*d(2)) / (e(1)/e(2)*w(2) - w(1));
        alpha = d(2)/e(2) + w(2)/e(2)*beta;
    end
    if and(alpha>=0,alpha<=1)  % intersection is within edge
        n_IntsP2 = n_IntsP2 + 1;
        IntsP2(1:2,n_IntsP2) = v1 + alpha*e - M;     % intersection point relative to M
        if beta >= 0
            r_IntsP2(1,n_IntsP2) = + norm(IntsP2(1:2,n_IntsP2),2); % radius of intersection point
        else
            r_IntsP2(1,n_IntsP2) = - norm(IntsP2(1:2,n_IntsP2),2); % radius of intersection point
        end
        if and(r_IntsP2(1,n_IntsP2)>=r(1),r_IntsP2(1,n_IntsP2)<=r(2))
            ints = true;
            coll = true;
        end
    end
    
    if ~ints  % no intersection with cone lines
        alpha = norm(e,2);
        e = e / alpha;
        beta = (M(1)-v1(1))*e(1) + (M(2)-v1(2))*e(2);
        proj = v1 - M + beta*e;    % projection of M onto the edge, relative to M
        r_proj = norm(proj,2);
        w(1) = min([alpha;max([0;beta])]);
        proj_m = v1 - M + w(1)*e;  % projection of M onto the edge, relative to M
        r_proj_m = norm(proj_m,2);
        
% 2.2) Projection point in [r(1),r(2)]
  
        if and( r_proj_m>=r(1) , r_proj_m<=r(2) ) % calculate psi of projection point 
            
        	if     proj_m(1) > 0
                psi = atan(proj_m(2)/proj_m(1));
            elseif proj_m(1) < 0
                psi = atan(proj_m(2)/proj_m(1)) + pi;
            else
                if     proj_m(2) > 0
                    psi = + pi/2;
                elseif proj_m(2) < 0
                    psi = - pi/2;
                else
                    psi = 0;
                end
            end
            psi = psi - floor((psi-phi(1))/2/pi)*2*pi;  % psi in [phi(1),phi(1)+2*pi)            
            if and( psi>=phi(1) , psi<=phi(2) ) % projection point is in the cone
                coll = true;
            end
        end
            
 % 2.3) Projection point in [0,r(1)]
        
        if and( r_proj<=r(1) , ~coll )
            
            if     proj(1) > 0
                psi = atan(proj(2)/proj(1));
            elseif proj(1) < 0
                psi = atan(proj(2)/proj(1)) + pi;
            else
                if     proj(2) > 0
                    psi = + pi/2;
                elseif proj(2) < 0
                    psi = - pi/2;
                else
                    psi = 0;
                end
            end
            r1 = norm(v1-M,2);  % radius of vertex 1
            r2 = norm(v2-M,2);  % radius of vertex 2
            
            % Vertex 1: check (r1 < r(1)) iff r_proj/r(1) in [w(1),w(2)]
            
            psi = psi - floor((psi-phi(1))/2/pi)*2*pi;  % psi in [phi(1),phi(1)+2*pi)
            if psi <= phi(1) + pi/2  % set lower bound
                w(1) = cos(psi-phi(1));
            elseif psi <= phi(2) + pi/2
                w(1) = 0;
            else
                w(1) = +1;
            end
            if psi <= phi(2)  % set upper bound
                w(2) = +1;
            elseif psi <= phi(2) + pi/2
                w(2) = cos(phi(2)-psi);
            else
                w(2) = -1;
            end
            if and( r_proj/r(1)>=w(1) , r_proj/r(1)<=w(2) )
                if and( beta>=0 , r1>=r(1) )
                    if or( beta<=alpha , r2<=r(2) )
                        ints = true;
                        coll = true;
                    end
                end
            end
            if psi <= phi(2)-3/2*pi  % additional upper bound
                if r_proj/r(1) <= cos(phi(2)-psi)
                    if and( beta>=0 , r1>=r(1) )
                        if or( beta<=alpha , r2<=r(2) )
                            ints = true;
                            coll = true;
                        end
                    end
                end
            end 

            % Vertex 2: check (r2 < r(1)) iff r_proj/r(1) in [w(1),w(2)]
            
            psi = psi - ceil((psi-phi(2))/2/pi)*2*pi;  % psi in (phi(2)-2*pi,phi(2)]
            if psi >= phi(1)   % set upper bound
                w(2) = +1;
            elseif psi >= phi(1)-pi/2
                w(2) = cos(psi-phi(1));
            else
                w(2) = -1;
            end
            if psi >= phi(2)-pi/2  % set lower bound
                w(1) = cos(psi-phi(2));
            elseif psi >= phi(1)-pi/2
                w(1) = 0;
            else
                w(1) = +1;
            end
            if and( r_proj/r(1)>=w(1) , r_proj/r(1)<=w(2) )
                if and( beta<=alpha , r2>=r(1) )
                    if or( beta>=0 , r1<=r(2) )
                        ints = true;
                        coll = true;
                    end
                end
            end
            if psi >= phi(1) + 3/2*pi  % additional upper bound
                if r_proj/r(1) <= cos(psi-phi(1))
                    if and( beta<=alpha , r2>=r(1) )
                        if or( beta>=0 , r1<=r(2) )
                            ints = true;
                            coll = true;
                        end
                    end
                end
            end 
            
        end
    end
    
    if coll
        break  % collision already detected
    end
end
 
% 2.4) Full inlusion by the obstacle
    
if ~coll  % check for full inclusion, and set coll = true 
    if and(n_IntsP1>=2,n_IntsP2>=2)
        count_in = false;
        count_out = false;
        if n_IntsP1 > 2             % debug only!!!!
            display(['Ahhh 1!'])
            n_IntsP1 = 2;
            term = true;
        end
        for k = 1:n_IntsP1
            if     r_IntsP1(1,k) <= r(1)
                count_in = true;
            elseif r_IntsP1(1,k) >= r(2)
                count_out = true;
            else  % debug only!!
                error('There should be no intersection in [r(1),r(2)]!!')  % debug only!!
            end
        end
        if and(count_in,count_out)
            count_in = false;
            count_out = false;
            if n_IntsP2 > 2             % debug only!!!!
                display(['Ahhh 1!'])
                n_IntsP2 = 2;
                term = true;
            end
            for k = 1:n_IntsP2
               if     r_IntsP2(1,k) <= r(1)
                   count_in = true;
               elseif r_IntsP2(1,k) >= r(2)
                   count_out = true;
               else  % debug only!!
                   error('There should be no intersection in [r(1),r(2)]!!')  % debug only!!
               end
            end
            if and(count_in,count_out)
                coll = true;
            end
        end
    end
end

% 4) Plot (Debugging) ------------------------------------------------------------------------------

if 0
    
    figure()
    hold on
    axis equal

% 3.1) Plot obstacle (green / red)

    plotopt.linestyle = '-';
    plotopt.linewidth = 1.5;
    if coll
        plotopt.color = [1 0 0]; % red
    else
        plotopt.color = [0 1 0]; % red
    end
    plotopt.marker = '.';
    plotopt.markersize = 10;
    plotopt.markeredgecolor = [0 0 0];
    plotopt.markerfacecolor = [0 0 0];
    plot([Obs(:,1);Obs(1,1)],[Obs(:,2);Obs(1,2)],plotopt);

% 3.2) Plot cone (black)

    plotopt.linestyle = '-';
    plotopt.linewidth = 1.5;
    plotopt.color = [0 0 0];
    plotopt.markeredgecolor = [0 0 0];
    plotopt.markerfacecolor = [0 0 0];
    plotopt.marker = '*';
    plotopt.markersize = 3;
    plot([M(1)],[M(2)],plotopt)
    plotopt.linestyle = '-';
    plotopt.linewidth = 1.5;
    plotopt.color = [0 0 0];
    plotopt.markeredgecolor = [0 0 0];
    plotopt.markerfacecolor = [0 0 0];
    plotopt.marker = 'none';
    plotopt.markersize = 3;
    plot([M(1)+r(1)*cos(phi(1)), M(1)+r(2)*cos(phi(1))],[M(2)+r(1)*sin(phi(1)), M(2)+r(2)*sin(phi(1))],plotopt);
    plot([M(1)+r(1)*cos(phi(2)), M(1)+r(2)*cos(phi(2))],[M(2)+r(1)*sin(phi(2)), M(2)+r(2)*sin(phi(2))],plotopt);

% 3.3) Plot circle (black)

    cosphi = cos(phi(1)+linspace(0,1,1001)*(phi(2)-phi(1)));
    sinphi = sin(phi(1)+linspace(0,1,1001)*(phi(2)-phi(1)));
    plot([M(1)+r(1)*cosphi],[M(2)+r(1)*sinphi],plotopt);
    plot([M(1)+r(2)*cosphi],[M(2)+r(2)*sinphi],plotopt);

% 3.4) Plot intersection points (blue)

    plotopt.linestyle = ':';
    plotopt.linewidth = 1.5;
    plotopt.color = [0 0 1];
    plotopt.markeredgecolor = [0 0 1];
    plotopt.markerfacecolor = [0 0 1];
    plotopt.marker = '*';
    plotopt.markersize = 3;
    if n_IntsP1 > 0
        plot(M(1)+[0,IntsP1(1,1:n_IntsP1)],M(2)+[0,IntsP1(2,1:n_IntsP1)],plotopt)
    end
    if n_IntsP2 > 0
        plot(M(1)+[0,IntsP2(1,1:n_IntsP2)],M(2)+[0,IntsP2(2,1:n_IntsP2)],plotopt)
    end

end

return


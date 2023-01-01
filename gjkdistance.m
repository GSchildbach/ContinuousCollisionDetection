% Georg Schildbach, 03/May/2021 --- RRT+CCD Algorithm
% Implementation of the Gilbert-Johnson-Keerthi (GJK) distance algorithm in R^2
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
% - V-polytope 1: vertices of X in R^2 (2 x nX)
% - V-polytope 2: vertices of Y in R^2 (2 x nY)
% --------------------------------------------------------------------------------------------------
% OUTPUTS:
% - flag: indicates an overlap between X and Y (true / false)
% - dist: distance between X and Y (= 0 if there is a collision)
% --------------------------------------------------------------------------------------------------

function [flag,dist] = gjkdistance(X,Y)

% 1) Inputs  ---------------------------------------------------------------------------------------

% 1.1) Parameters

tol = 100*eps;  % eps: machine epsilon
maxit = 100;    % maximum number of iterations

% 1.2) Initialize

flag = true;    % collision flag
dist = 0;       % distance between X and Y
nX = uint8(0);
nY = uint8(0);
v = zeros(2,1);
w = zeros(2,1);
lambda = zeros(3,1);
D = zeros(3,1);

% 2) GJK-Algorithm ---------------------------------------------------------------------------------

% 2.1) Initialize 

nS = uint8(1);                        % number of nodes in the simplex
S = [X(1:2,1)-Y(1:2,1) , zeros(2,2)]; % first nodes

for k = 1:maxit
    
    %display(['*** Iteration number ', num2str(k,'%d') ' ***'])
    
% 2.2) Simplex consisting of one node
    
    if     nS == 1  
        if and( abs(S(1,1))<=tol , abs(S(2,1))<=tol )
            flag = true;
            dist = 0;
            break  % optimum found and equal to zero
        else
            nX = suppfct(-S(1:2,1),X);
            nY = suppfct(S(1:2,1),Y);
            v = X(1:2,nX) - Y(1:2,nY);
            if and( abs(S(1,1)-v(1,1))<=tol , abs(S(2,1)-v(2,1))<=tol )
                flag = false;
                dist = norm(v(1:2),2);
                break  % optimum found and not equal to zero
            else
                nS = 2;
                S(1:2,2) = v; % second node
            end
        end
        
% 2.3) Simplex consisting of two nodes
 
    elseif nS == 2
        w = S(1:2,2) - S(1:2,1);
        lambda(2,1) = norm(w,2);
        if lambda(2,1) == 0
            nS = 1;
        else
            w = w / lambda(2,1);
            lambda(1,1) = - S(1:2,1)' * w;
            if     lambda(1,1) < 0
                lambda(1,1) = 0;
            elseif lambda(1,1) > lambda(2,1)
                lambda(1,1) = lambda(2,1);
            end
            w = S(1:2,1) + lambda(1,1)*w;  % closest point on the node connection to the origin
            if and( abs(w(1,1))<=tol , abs(w(2,1))<=tol )
                flag = true;
                dist = 0;
                break  % optimum found and equal to zero
            else
                nX = suppfct(-w(1:2),X);
                nY = suppfct(w(1:2),Y);
                v = X(1:2,nX) - Y(1:2,nY);
                if     and( abs(v(1,1)-S(1,1))<=tol , abs(v(2,1)-S(2,1))<=tol )
                    flag = false;
                    dist = norm(w(1:2),2);
                    break  % optimum found and not equal to zero
                elseif and( abs(v(1,1)-S(1,2))<=tol , abs(v(2,1)-S(2,2))<=tol )
                    flag = false;
                    dist = norm(w(1:2),2);
                    break  % optimum found and not equal to zero
                else
                    nS = 3;
                    S(1:2,3) = v; % third node
                end
            end
        end
        
% 2.4) Simplex consisting of three nodes

    elseif nS == 3
        D(1,1) = S(1,1)*S(2,2) - S(1,2)*S(2,1);
        D(2,1) = S(1,1)*S(2,3) - S(1,3)*S(2,1);
        D(3,1) = S(1,2)*S(2,3) - S(1,3)*S(2,2);
        if     D(1,1) ~= 0
            lambda(1,1) = - D(3,1) / D(1,1);
            lambda(2,1) = D(2,1) / D(1,1);
            lambda(3,1) = -1;
        elseif D(2,1) ~= 0
            lambda(1,1) = D(3,1) / D(2,1);
            lambda(2,1) = -1;
            lambda(3,1) = D(1,1) / D(2,1);
        elseif D(3,1) ~= 0
            lambda(1,1) = -1;
            lambda(2,1) = D(2,1) / D(3,1);
            lambda(3,1) = - D(1,1) / D(3,1);
        else  % points are co-linear
            if     or( S(1,1)==S(1,3) , S(2,1)==S(2,3) )
                nS = 2;
            elseif or( S(1,2)==S(1,3) , S(2,2)==S(2,3) )
                nS = 2;
            elseif or( S(1,1)==S(1,2) , S(2,1)==S(2,2) )
                nS = 2;
                S(1:2,1) = S(1:2,3);
            else
                lambda(1,1) = (S(1,1)-S(1,2)) / (S(1,3)-S(1,2));
                lambda(2,1) = (S(1,2)-S(1,3)) / (S(1,1)-S(1,3));
                lambda(3,1) = (S(1,3)-S(1,1)) / (S(1,2)-S(1,1));
                if     and( lambda(1,1)>=0 , lambda(1,1)<=1 )
                    nS = 2;
                    S(1:2,1) = S(1:2,3);
                elseif and( lambda(2,1)>=0 , lambda(2,1)<=1 )
                    nS = 2;
                    S(1:2,2) = S(1:2,3);
                elseif and( lambda(3,1)>=0 , lambda(3,1)<=1 )
                    nS = 2;
                end
            end
        end
        if nS == 3
            D(1,1) = sum(lambda(:,1));
            lambda(1,1) = lambda(1,1) / D(1,1);
            lambda(2,1) = lambda(2,1) / D(1,1);
            lambda(3,1) = lambda(3,1) / D(1,1);
            if and( lambda(1,1)>-tol , and( lambda(2,1)>-tol , lambda(3,1)>-tol ))
                flag = true;
                dist = 0;
                break  % optimum found and equal to zero
            else
                v = S(1:2,2) - S(1:2,1);
                lambda(2,1) = norm(v,2);
                v = v / lambda(2,1);
                lambda(1,1) = - S(1:2,1)' * v;
                if     lambda(1,1) < 0
                    lambda(1,1) = 0;
                elseif lambda(1,1) > lambda(2,1)
                    lambda(1,1) = lambda(2,1);
                end
                v = S(1:2,1) + lambda(1,1)*v;  % closest point on the node connection to the origin
                lambda(3,1) = norm(v,2);
                w = v;
                v = S(1:2,3) - S(1:2,2);
                lambda(2,1) = norm(v,2);
                v = v / lambda(2,1);
                lambda(1,1) = - S(1:2,2)' * v;
                if     lambda(1,1) < 0
                    lambda(1,1) = 0;
                elseif lambda(1,1) > lambda(2,1)
                    lambda(1,1) = lambda(2,1);
                end
                v = S(1:2,2) + lambda(1,1)*v;  % closest point on the node connection to the origin
                lambda(2,1) = norm(v,2);
                if lambda(2,1) < lambda(3,1)
                    lambda(3,1) = lambda(2,1);
                    w = v;
                    nS = 1;  % farthest node
                end
                v = S(1:2,1) - S(1:2,3);
                lambda(2,1) = norm(v,2);
                v = v / lambda(2,1);
                lambda(1,1) = - S(1:2,3)' * v;
                if     lambda(1,1) < 0
                    lambda(1,1) = 0;
                elseif lambda(1,1) > lambda(2,1)
                    lambda(1,1) = lambda(2,1);
                end
                v = S(1:2,3) + lambda(1,1)*v;  % closest point on the node connection to the origin
                lambda(2,1) = norm(v,2);
                if lambda(2,1) < lambda(3,1)
                    lambda(3,1) = lambda(2,1);
                    w = v;
                    nS = 2;  % farthest node
                end
                nX = suppfct(-w(1:2),X);
                nY = suppfct(w(1:2),Y);
                v = X(1:2,nX) - Y(1:2,nY);
                if     and( abs(v(1,1)-S(1,1))<=tol , abs(v(2,1)-S(2,1))<=tol )
                    flag = false;
                    dist = norm(w(1:2),2);
                    break  % optimum found and not equal to zero
                elseif and( abs(v(1,1)-S(1,2))<=tol , abs(v(2,1)-S(2,2))<=tol )
                    flag = false;
                    dist = norm(w(1:2),2);
                    break  % optimum found and not equal to zero
                elseif and( abs(v(1,1)-S(1,3))<=tol , abs(v(2,1)-S(2,3))<=tol )
                    flag = false;
                    dist = norm(w(1:2),2);
                    break  % optimum found and not equal to zero
                else
                    S(1:2,nS) = v; % replace farthest node
                    nS = 3;
                end
            end
        end
    end
    
end

return


% ==================================================================================================


% --------------------------------------------------------------------------------------------------
% Support function: vertex = suppfct(v,X)
% --------------------------------------------------------------------------------------------------
% INPUTS:
% - V-polytope 1: vertices of X in R^2 (2 x nX)
% - vector v (2 x 1)
% --------------------------------------------------------------------------------------------------
% OUTPUTS:
% - vertex k that maximizes the support function argmax(v'*X(:,k))
% --------------------------------------------------------------------------------------------------

function vertex = suppfct(v,X)

supp = -Inf;
vertex = uint8(0);

for k = 1:1:size(X,2)
    check = v' * X(1:2,k);
    if check > supp
        supp = check;
        vertex = k;
    end
end

return





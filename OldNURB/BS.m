function [N, M] = BS(p1, p2, n1,n2, t1, t2, xis, etas)
% Given knot vectors, control points, weights, and points to evaluate (in parameter space),
% order of functions, evaluates the basis functions at xi/eta
% To be general, xi here is just a stencil from 0 -> 1 that contains
% fractions on which to evaluate xi along the whole length of the knot
% vector

%% Evaluate the 1st set of basis functions

% iter is number of basis functions to evaluate for each N
iter = n1+p1;
% N contains all basis functions of all orders
N = cell(iter,p1);
% Evaluation points
ev1 = length(xis);
% tmp is the temporary holder for the evaluation of a given basis function
tmp = zeros(ev1,1);
% basis functions of 1 order lower at i and i+1
tmpm=tmp;
tmpp=tmp;

% Convert stencil to actual xi's
low = t1(p1+1);
high= t1(n1+1);
xi = (high - low)*xis + low;% CHANGE BACK FOR COMP???

% Outer loop goes from order 0 -> order p basis functions
for i = 1:p1+1
%   Inner loop evaluates the basis functions
    for j = 1:iter
%       Get lower order basis functions
        if i>1
            tmpm = N{j,i-1};
            tmpp = N{j+1,i-1};
        end
%       Loop over knot values to eval
        for kk = 1:ev1
%               Evaluate 0th order basis fns
                if i==1
                    if xi(kk)< t1(j+1)&& xi(kk)>=t1(j)
                        tmp(kk) = 1;
                    else
                        tmp(kk) = 0;
                    end
%                   Deal with the very last point to avoid NaNs
                    if xi(kk)==t1(end) && tmp(kk-1,ll) == 1
                        tmp(kk)=1;
                    end
%               Evaluate the ith order basis function    
                else
                    tmp(kk) = 0;
%                   This is to deal with repeated knots
                    if t1(j+i-1)~=t1(j)
                        tmp(kk) = (xi(kk) - t1(j)) ... 
                                      /(t1(j+i-1) - t1(j)  )*tmpm(kk);
                    end
                    if t1(j+i) ~= t1(j+1)
                        tmp(kk) = tmp(kk) + (t1(j+i) - xi(kk)) ...
                                     /(t1(j+i  ) - t1(j+1))*tmpp(kk);
                    end
                end
        end
        N{j,i} = tmp;
    end
    iter = iter-1;
end

%% Now for the 2nd set of basis fns

% iter is number of basis functions to evaluate for each M
iter = n2+p2;
% N contains all basis functions of all orders
M = cell(iter,p2);
% tmp is the temporary holder for the evaluation of a given basis function
ev2 = length(etas);
tmp = zeros(ev2,1);
% basis functions of 1 order lower at i and i+1
tmpm=tmp;
tmpp=tmp;

% Convert stencil to actual eta's
low = t2(p2+1);
high= t2(n2+2);
eta = (high - low)*etas + low; %CHANGE BACK AGAIN!

% Outer loop goes from order 0 -> order P basis functions
for i = 1:p2+1
%   Inner loop evaluates the basis functions
    for j = 1:iter
%       Get lower order basis functions
        if i>1
            tmpm = M{j,i-1};
            tmpp = M{j+1,i-1};
        end
%       Loop over knot values to eval
        for kk = 1:ev2(1)
%               Evaluate 0th order basis fns
                if i==1
                    if eta(kk)< t2(j+1)&& eta(kk)>=t2(j)
                        tmp(kk) = 1;
                    else
                        tmp(kk) = 0;
                    end
%                   Deal with the very last point to avoid NaNs
                    if eta(kk)==t2(end) && tmp(kk-1,ll) == 1
                        tmp(kk)=1;
                    end
%               Evaluate the ith order basis function    
                else
                    tmp(kk) = 0;
%                   This is to deal with repeated knots
                    if t2(j+i-1)~=t2(j)
                        tmp(kk) = (eta(kk) - t2(j)  ) ...
                                     /(t2(j+i-1) - t2(j)  )*tmpm(kk);
                    end
                    if  t2(j+i) ~= t2(j+1)
                        tmp(kk) = tmp(kk) + (t2(j+i) - eta(kk))...
                                     /(t2(j+i  ) - t2(j+1))*tmpp(kk);
                    end
                end
        end
        M{j,i} = tmp;
    end
    iter = iter-1;
end

end
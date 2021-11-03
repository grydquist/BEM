function dN = BSd(N, t, xis)
% Given the big NURBS cell arrays, calculates the derivatives at
% those same points,


%% Compute the derivatives just with respect to the basis function here
% info about the curves
[iter, p] = size(N);
p = p - 1;
n = iter - p;

% resizing xi
low = t(p+1);
high= t(n+1);
xi = xis;%(high - low)*xis + low; CHANGE BACK FOR OTHER
dN = cell(n,1);

ev = size(xi);
tmp = zeros(ev(1),ev(2));

% Loop over basis functions
for i = 1:n
%   Evaluate the derivative of the ith basis function at xi(j)
    Np1  = N{i  ,p};
    Ni1p1= N{i+1,p};
    
%   Actually calculate derivatives at xi (2.12)
    tmp(:) = 0;
    for j1 = 1:ev(1)
        for j2 = 1:ev(2)
            if t(i+p) ~= t(i)
                tmp(j1,j2) = p/(t(i+p) - t(i))*Np1(j1,j2);
            end
            if t(i+p+1) ~= t(i+1)
                tmp(j1,j2) = tmp(j1,j2) - p/(t(i+p+1)-t(i+1))*Ni1p1(j1,j2);
            end
        end
    end
    dN{i} = tmp;
end 

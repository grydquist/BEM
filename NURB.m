function R = NURB(N,M,w)
% Given basis functions evaluated at known xi and eta as well as control
% points, evaluates X, locations of points at (xi,eta) on surface
% The weights are the fourth row of P

%% problem info
ev1 = size(N{1,1});
ev2 = size(M{1,1});

[iter,p1] = size(N);
p1 = p1 - 1;
n1 = iter-p1;

[iter,p2] = size(M);
p2 = p2-1;
n2 = iter-p2;

% actual NURBS storage
R = cell(n1,n2);
tmp5 = zeros(ev1(1),ev2(1),ev1(2),ev2(2));

%% Weighted denominator

% Just getting weighted denominator
% Size is Nel1 x Nel2 x ng1 x ng2
NMw = zeros(ev1(1),ev2(1),ev1(2),ev2(2));
for i = 1:n1
    tmp1 = N{i,end};
    for j = 1:n2
        tmp2 = M{j,end};
        for k1 = 1:ev1(1)
            for k2 = 1:ev1(2)
                for l1 = 1:ev2(1)
                    for l2 = 1:ev2(2)
                        NMw(k1,l1,k2,l2) = NMw(k1,l1,k2,l2) ...
                        + tmp1(k1,k2)*tmp2(l1,l2)*w(i,j);
                    end
                end
            end
        end
    end
end

%% Actual evaluation of NURBS at locations
% Loop over xi basis functions at all xi
for i = 1:n1
    tmp1 = N{i,end};
%   Loop over basis functions at all eta
    for j = 1:n2
        tmp2 = M{j,end};
%       Loop over xis at a element/GP combo
        for k1 = 1:ev1(1)
            for k2 = 1:ev1(2)
%               Loop over etas at a element/GP combo
                for l1 = 1:ev2(1)
                    for l2 = 1:ev2(2)
%                       Evaluate NURBS basis here
                        tmp5(k1,l1,k2,l2) = ...
                         tmp1(k1,k2)*tmp2(l1,l2)*w(i,j)/NMw(k1,l1,k2,l2);
                    end
                end
            end
        end
%       Put back into NURBS basis cell
        R{i,j} = tmp5;
    end
end

end 
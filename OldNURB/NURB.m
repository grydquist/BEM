function X = NURB(N,M,P)
% Given basis functions evaluated at known xi and eta as well as control
% points, evaluates X, locations of points at (xi,eta) on surface
% The weights are the fourth row of P

%% problem info
ev1 = length(N{1,1});
ev2 = length(M{1,1});

[iter,p1] = size(N);
p1 = p1 - 1;
n1 = iter-p1;

[iter,p2] = size(M);
p2 = p2-1;
n2 = iter-p2;
%% Weighted denominator

% Size is Nel1 x Nel2 x ng1 x ng2
NMw = zeros(ev1,ev2);
for i = 1:n1
    tmp1 = N{i,end};
    for j = 1:n2
        tmp2 = M{j,end};
        for k1 = 1:ev1
            for l1 = 1:ev2
                NMw(k1,l1) = NMw(k1,l1) ...
                + tmp1(k1)*tmp2(l1)*P(4,i,j);
            end
        end
    end
end

%% Evaluate at a bunch of points
% Size is Nel1 x Nel2 x ng1 x ng2
uu = zeros(ev1,ev2);
vv = uu;
ww = vv;

% Loop over first set of basis functions
for i =1:n1
%   ith basis function of N
    tmp1 = N{i,end};
    for j = 1:n2
%       jth basis function of M
        tmp2 = M{j,end};
%       Loop over evaluations on eta and xi
        for k1 = 1:ev1
            for l1 = 1:ev2
                uu(k1,l1) = uu(k1,l1) ...
                   + tmp1(k1)*tmp2(l1)*P(1,i,j)*P(4,i,j)...
                   /NMw(k1,l1);
                vv(k1,l1) = vv(k1,l1) ...
                   + tmp1(k1)*tmp2(l1)*P(2,i,j)*P(4,i,j)...
                   /NMw(k1,l1);
                ww(k1,l1) = ww(k1,l1) ...
                   + tmp1(k1)*tmp2(l1)*P(3,i,j)*P(4,i,j)...
                   /NMw(k1,l1);
            end            
        end
    end
end

% Size is nsd x Nel1 x Nel2 x ng1 x ng2
X = zeros(3,ev1,ev2);
X(1,:,:) = uu;
X(2,:,:) = vv;
X(3,:,:) = ww;

end 
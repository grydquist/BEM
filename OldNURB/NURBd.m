function [dXxi, dXeta] = NURBd(N,M,dN,dM,P)
% Evaluates derivatives on surface with respect to both eta and xi

%% problem info
ev1 = length(N{1,1});
ev2 = length(M{1,1});

[iter,p1] = size(N);
p1 = p1 - 1;
n1 = iter-p1;

[iter,p2] = size(M);
p2 = p2-1;
n2 = iter-p2;

%% Weighted denominators, including those with derivatives

NMw = zeros(ev1,ev2);
dNMw= NMw;
NdMw= NMw;
for i = 1:n1
    tmp1 = N{i,end};
    tmp3 =dN{i};
    for j = 1:n2
        tmp2 = M{j,end};
        tmp4 =dM{j};
        for k = 1:ev1
            for l = 1:ev2
                NMw(k,l) = NMw(k,l) ...
                    + tmp1(k)*tmp2(l)*P(4,i,j);
                dNMw(k,l)= dNMw(k,l)...
                    + tmp3(k)*tmp2(l)*P(4,i,j);
                NdMw(k,l)= NdMw(k,l)...
                    + tmp1(k)*tmp4(l)*P(4,i,j);
            end
        end
    end
end

%% Evaluate at a bunch of points
% Size is Nel1 x Nel2 x ng1 x ng2
uuxi = zeros(ev1,ev2);
vvxi = uuxi;
wwxi = vvxi;

uueta = zeros(ev1,ev2);
vveta = uueta;
wweta = vveta;

% Loop over first set of basis functions
for i =1:n1
%   ith basis function of N
    tmp1 = N{i,end};
    tmp3 =dN{i,end};
    for j = 1:n2
%       jth basis function of M
        tmp2 = M{j,end};
        tmp4 =dM{j,end};
%       Loop over evaluations on eta and xi
        for k = 1:ev1
            for l = 1:ev2
                multxi = (tmp3(k)*tmp2(l)*NMw(k,l)...
                    - tmp1(k)*tmp2(l)*dNMw(k,l)) ...
                    / (NMw(k,l)^2)*P(4,i,j);
                uuxi(k,l) = uuxi(k,l) + multxi * P(1,i,j);
                vvxi(k,l) = vvxi(k,l) + multxi * P(2,i,j);
                wwxi(k,l) = wwxi(k,l) + multxi * P(3,i,j);

                multeta= (tmp1(k)*tmp4(l)*NMw(k,l)...
                    - tmp1(k)*tmp2(l)*NdMw(k,l)) ...
                    / (NMw(k,l)^2)*P(4,i,j);
                uueta(k,l) = uueta(k,l) + multeta * P(1,i,j);
                vveta(k,l) = vveta(k,l) + multeta * P(2,i,j);
                wweta(k,l) = wweta(k,l) + multeta * P(3,i,j);
            end
        end
    end
end

%% Plot points
% reshape into vectors
dXeta = zeros(3,ev1,ev2);
dXeta(1,:,:) = uueta;
dXeta(2,:,:) = vveta;
dXeta(3,:,:) = wweta;

dXxi = zeros(3,ev1,ev2);
dXxi(1,:,:) = uuxi;
dXxi(2,:,:) = vvxi;
dXxi(3,:,:) = wwxi;


end 
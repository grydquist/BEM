function [dXxi, dXeta] = NURBd(N,M,dN,dM,P,w)
% Evaluates derivatives on surface with respect to both eta and xi

%% problem info
ev1 = size(N{1,1});
ev2 = size(M{1,1});

[iter,p1] = size(N);
p1 = p1 - 1;
n1 = iter-p1;

[iter,p2] = size(M);
p2 = p2-1;
n2 = iter-p2;

%% Weighted denominators, including those with derivatives

% Size is Nel1 x Nel2 x ng1 x ng2
NMw = zeros(ev1(1),ev2(1),ev1(2),ev2(2));
dNMw= NMw;
NdMw= NMw;
for i = 1:n1
    tmp1 = N{i,end};
    tmp3 =dN{i};
    for j = 1:n2
        tmp2 = M{j,end};
        tmp4 =dM{j};
        for k1 = 1:ev1(1)
            for k2 = 1:ev1(2)
                for l1 = 1:ev2(1)
                    for l2 = 1:ev2(2)
                        NMw(k1,l1,k2,l2) = NMw(k1,l1,k2,l2) ...
                            + tmp1(k1,k2)*tmp2(l1,l2)*w(i,j);
                        dNMw(k1,l1,k2,l2)= dNMw(k1,l1,k2,l2)...
                            + tmp3(k1,k2)*tmp2(l1,l2)*w(i,j);
                        NdMw(k1,l1,k2,l2)= NdMw(k1,l1,k2,l2)...
                            + tmp1(k1,k2)*tmp4(l1,l2)*w(i,j);
                    end
                end
            end
        end
    end
end

%% Evaluate at a bunch of points
% Size is Nel1 x Nel2 x ng1 x ng2
uuxi = zeros(ev1(1),ev2(1),ev1(2),ev2(2));
vvxi = uuxi;
wwxi = vvxi;

uueta = zeros(ev1(1),ev2(1),ev1(2),ev2(2));
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
        for k1 = 1:ev1(1)
            for k2 = 1:ev1(2)
                for l1 = 1:ev2(1)
                    for l2 = 1:ev2(2)
                        multxi = (tmp3(k1,k2)*tmp2(l1,l2)*NMw(k1,l1,k2,l2)...
                            - tmp1(k1,k2)*tmp2(l1,l2)*dNMw(k1,l1,k2,l2)) ...
                            / (NMw(k1,l1,k2,l2)^2)*w(i,j);
                        uuxi(k1,l1,k2,l2) = uuxi(k1,l1,k2,l2) + multxi * P(1,i,j);
                        vvxi(k1,l1,k2,l2) = vvxi(k1,l1,k2,l2) + multxi * P(2,i,j);
                        wwxi(k1,l1,k2,l2) = wwxi(k1,l1,k2,l2) + multxi * P(3,i,j);

                        multeta= (tmp1(k1,k2)*tmp4(l1,l2)*NMw(k1,l1,k2,l2)...
                            - tmp1(k1,k2)*tmp2(l1,l2)*NdMw(k1,l1,k2,l2)) ...
                            / (NMw(k1,l1,k2,l2)^2)*w(i,j);
                        uueta(k1,l1,k2,l2) = uueta(k1,l1,k2,l2) + multeta * P(1,i,j);
                        vveta(k1,l1,k2,l2) = vveta(k1,l1,k2,l2) + multeta * P(2,i,j);
                        wweta(k1,l1,k2,l2) = wweta(k1,l1,k2,l2) + multeta * P(3,i,j);
                    end
                end
            end            
        end
    end
end

%% Plot points
% reshape into vectors
dXeta = zeros(3,ev1(1),ev2(1),ev1(2),ev2(2));
dXeta(1,:,:,:,:) = uueta;
dXeta(2,:,:,:,:) = vveta;
dXeta(3,:,:,:,:) = wweta;

dXxi = zeros(3,ev1(1),ev2(1),ev1(2),ev2(2));
dXxi(1,:,:,:,:) = uuxi;
dXxi(2,:,:,:,:) = vvxi;
dXxi(3,:,:,:,:) = wwxi;


end 
function X = NURBev(R,P)
% Given a set of NURBS basis functions and control points, evaluates the
% location of those control points

%% problem info
et = size(R{1});
if length(et) == 2
    et(3) = 1;
    et(4) = 1;
end

if length(et) == 3
    et(4) = 1;
end

ev1 = [et(1),et(3)];
ev2 = [et(2),et(4)];

[n1,n2] = size(R);

%% Evaluate at a bunch of points
% Size is Nel1 x Nel2 x ng1 x ng2
uu = zeros(ev1(1),ev2(1),ev1(2),ev2(2));
vv = uu;
ww = vv;

% Big loop over basis functions
for i =1:n1
    for j = 1:n2
%       Current basis function evaluated at all xi/eta
        tmp5 = R{i,j};
%       Element/Gauss point loop over xi
        for k1 = 1:ev1(1)
            for k2 = 1:ev1(2)
%               Element/Gauss point loop over eta
                for l1 = 1:ev2(1)
                    for l2 = 1:ev2(2)
                        uu(k1,l1,k2,l2) = uu(k1,l1,k2,l2) + ...
                            tmp5(k1,l1,k2,l2)*P(1,i,j);
                        vv(k1,l1,k2,l2) = vv(k1,l1,k2,l2) + ...
                            tmp5(k1,l1,k2,l2)*P(2,i,j);
                        ww(k1,l1,k2,l2) = ww(k1,l1,k2,l2) + ...
                            tmp5(k1,l1,k2,l2)*P(3,i,j);
                    end
                end
            end            
        end
    end
end

% Size is nsd x Nel1 x Nel2 x ng1 x ng2
X = zeros(3,ev1(1),ev2(1),ev1(2),ev2(2));
X(1,:,:,:,:) = uu;
X(2,:,:,:,:) = vv;
X(3,:,:,:,:) = ww;

end
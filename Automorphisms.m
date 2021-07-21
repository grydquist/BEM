L = eye(2);
clf
% plot([0,L(1,1)],[0,L(1,2)]);hold on;plot([0,L(2,1)],[0,L(2,2)])
M3 = [1,2;0,1];
M2 = [1,0;2,1];
M6 = [1,0;2,1];
M4 = [1,1;1,2];
M5 = [2,-2;-2,3];
M7 = [-1,2;0,1];
M8 = [1,1;0,1];
M9 = [1,0;1,1];

A = [1,1;1,-1];
L2 = L*expm(A*.65);
plot([0,L2(1,1)],[0,L2(1,2)],'--');hold on;plot([0,L2(2,1)],[0,L2(2,2)],'--')
plot([0,L2(1,1)]+ L2(2,1),[0,L2(1,2)]+L2(2,2),'--');hold on;plot([0,L2(2,1)]+L2(1,1),[0,L2(2,2)]+L2(1,2),'--')
% L3 = inv(M8)*L2;
% L3 = (M7)*L2;
L3 = inv(M3)*L2;
plot([0,L3(1,1)],[0,L3(1,2)]);hold on;plot([0,L3(2,1)],[0,L3(2,2)])
plot([0,L3(1,1)] + L3(2,1),[0,L3(1,2)] + L3(2,2) );hold on;plot([0,L3(2,1)]+L3(1,1),[0,L3(2,2)]+L3(1,2))

axis([-2,4,-2,2])
pbaspect([1.5,1,1])

% disp([abs(acos(dot(L2(:,1),L2(:,2))/(norm(L2(:,1))*norm(L2(:,2))))/pi*180-90), ...
%     abs(acos(dot(L3(:,1),L3(:,2))/(norm(L3(:,1))*norm(L3(:,2))))/pi*180-90)])

p = zeros(2,16);
p2= p;
p3= p;

p(:,1) = [1,.5] - 0*[L2(1,1),L2(2,1)] - 2*[L2(1,2),L2(2,2)];
p2(:,1) = [0.25,.75] - 0*[L2(1,1),L2(2,1)] - 2*[L2(1,2),L2(2,2)];
p3(:,1) = [0,0] - 0*[L2(1,1),L2(2,1)] - 2*[L2(1,2),L2(2,2)];

it = 1;
for i = 0:4
    for j = 0:4
        it = it + 1;
        p(:,it) = p(:,1)  + i*[L2(1,1);L2(2,1)];
        p(:,it) = p(:,it) + j*[L2(1,2);L2(2,2)];
        p2(:,it) = p2(:,1)  + i*[L2(1,1);L2(2,1)];
        p2(:,it) = p2(:,it) + j*[L2(1,2);L2(2,2)];
        p3(:,it) = p3(:,1)  + i*[L2(1,1);L2(2,1)];
        p3(:,it) = p3(:,it) + j*[L2(1,2);L2(2,2)];
        
    end
end

scatter(p(1,:),p(2,:))
scatter(p2(1,:),p2(2,:))
scatter(p3(1,:),p3(2,:))

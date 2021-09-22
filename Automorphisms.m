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
M10= [1,-2;-2,5];

A = [1,1;1,-1];
L2 = L*expm(A*1.2);
plot([0,L2(1,1)],[0,L2(1,2)],'--');hold on;plot([0,L2(2,1)],[0,L2(2,2)],'--')
plot([0,L2(1,1)]+ L2(2,1),[0,L2(1,2)]+L2(2,2),'--');hold on;plot([0,L2(2,1)]+L2(1,1),[0,L2(2,2)]+L2(1,2),'--')
% L3 = inv(M8)*L2;
% L3 = (M7)*L2;
L3 = inv(M3)*L2;
% L3 = M10*L2;
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
for i = -4:4
    for j = -4:4
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

%% Try by just using nearest lattice point
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'periodic.gif';

% This really seems to work?????
L2 = L;
LL = L;
mdot = [0.35;0.35];
mdot2 = mdot.*1.5;
% Move the lattice a little each time
maxt = 100;
theta = linspace(1,20,maxt);
uu = zeros(2,2,maxt);
uu(1,1,:) = sin(theta);
uu(1,2,:) = 10;
uu(2,1,:) = 8*cos(2*theta);
uu(2,2,:) = -sin(theta);
for t = 1:maxt
dt = 0.01;
A = uu(:,:,t);

L2 = L2 + A*L2*dt;
LL = LL + A*LL*dt;
mdot = mdot + A*mdot*dt;
mdot2 = mdot2 + A*mdot2*dt;

d = max([norm(L2(:,1)),norm(L2(:,2))]);
        if t==20
            disp(1)
        end
% Starting at the lattice point at 0,0, find the nearest lattice point
for i = 1:1
    for j = -1:1
        if i==0 && j==0; continue; end
        pp = i*[L2(1,1);L2(2,1)] + j*[L2(1,2);L2(2,2)];
        if d>=norm(pp)
            d = norm(pp);
            ii = i;
            jj = j;
        end
    end
end
% ii = 1;jj=0;
myp = ii*[L2(1,1);L2(2,1)] + jj*[L2(1,2);L2(2,2)];
L2(:,1) = myp;
% Then find the second nearest, as well as angle (don't want too sheared)
ang = pi/2;
for i = -1:1
    for j = 1:1
        if i==0 && j==0; continue; end
        pp = i*[L2(1,1);L2(2,1)] + j*[L2(1,2);L2(2,2)];
        ang2 = abs(pi/2 - acos(dot(pp,myp)/norm(pp)/norm(myp)));
        if ang>=ang2
            ang = ang2;
            ii2 = i;
            jj2 = j;
        end
    end
end
% ii2 = 0;jj2=1;
myp2 = ii2*[L2(1,1);L2(2,1)] + jj2*[L2(1,2);L2(2,2)];
L2 = [myp,myp2];
clf
L3 = L2;%[myp';myp2'];
plot([0,L3(1,1)],[0,L3(2,1)],'r');hold on;plot([0,L3(1,2)],[0,L3(2,2)],'r')
plot([0,L3(1,1)] + L3(1,2),[0,L3(2,1)] + L3(2,2),'r');hold on;plot([0,L3(1,2)]+L3(1,1),[0,L3(2,2)]+L3(2,1),'r')

plot([0,LL(1,1)],[0,LL(2,1)],'b');hold on;plot([0,LL(1,2)],[0,LL(2,2)],'b')
plot([0,LL(1,1)] + LL(1,2),[0,LL(2,1)] + LL(2,2),'b');hold on;plot([0,LL(1,2)]+LL(1,1),[0,LL(2,2)]+LL(2,1),'b')


for i = -5:5
    for j = -5:5
        it = it + 1;
        plot(mdot(1) + i*LL(1,1) + j*LL(1,2),mdot(2) + i*LL(2,1) + j*LL(2,2),'ko');
        plot(mdot2(1) + i*LL(1,1) + j*LL(1,2),mdot2(2) + i*LL(2,1) + j*LL(2,2),'go');
        
    end
end
% plot(mdot(1),mdot(2),'o');
% plot(mdot2(1),mdot2(2),'o');

axis([-2,2,-2,2])
drawnow
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if t == 1 ],'b')
% plot([0,LL(1,1)] + LL(1,2),[0,LL(2,1)] + LL(2,2),'b');hold on;plot([0,
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
% pause(0.1);
end
L3 = L2;%[myp';myp2'];
% plot([0,L3(1,1)],[0,L3(2,1)]);hold on;plot([0,L3(1,2)],[0,L3(2,2)])
% plot([0,L3(1,1)] + L3(1,2),[0,L3(2,1)] + L3(2,2) );hold on;plot([0,L3(1,2)]+L3(1,1),[0,L3(2,2)]+L3(2,1))
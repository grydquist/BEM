% Simpson et al was a massive help here
mu = 1;
U = [0;1;0];

%% Actual code - NURBS Construction
% Necessary NURBS info

% order
p1 = 2;
p2 = 2;

% number of control points/basis functions
n1 = 9;
n2 = 5;

% size of knot vector
k1 = n1+p1+1;
k2 = n2+p2+1;

% actual knot vector
t1 = [0,0,0,.25,.25,.5,.5,.75,.75,1,1,1]*4;
t2 = [0,0,0,0.5,0.5,1,1,1]*2;

% Control net

P = zeros(3,n1,n2);
w = zeros(n1,n2);
xs = [1,1,0,-1,-1];
ys  = [0,1,1,1,0];
r2 = sqrt(2)/2;
w1 = [1,r2,1,r2,1,r2,1,r2,1];
w2 = [1,r2,1,r2,1];

for i = 1:n1
    P(2,i,:) = xs;
    w(i,:) = w2*w1(i);
end

P(1,1,:) = ys;
P(1,2,:) = ys;
P(1,3,:) = 0;
P(1,4,:) = -ys;
P(1,5,:) = -ys;
P(1,6,:) = -ys;
P(1,7,:) = 0;
P(1,8,:) = ys;
P(1,9,:) = ys;

P(3,1,:) = 0;
P(3,2,:) = ys;
P(3,3,:) = ys;
P(3,4,:) = ys;
P(3,5,:) = 0;
P(3,6,:) = -ys;
P(3,7,:) = -ys;
P(3,8,:) = -ys;
P(3,9,:) = 0;

w(:,1) = 1;
w(:,5) = 1;

%% Pre-processing calculations
% Number of elements in each direction: this will change with periodic
% splines b/c we won't use all the knot spans, but for now this is just all
% the unique spans

% Get knot length/size of elements in parent domain
i = 1;
while(t1(i)==0)
    i=i+1;
end
d1 = t1(i);

i = 1;
while(t2(i)==0)
    i=i+1;
end
d2 = t2(i);

Nel1 = max(t1)/d1;
Nel2 = max(t2)/d2;

% Useful to have an ascending list without repeats (effectively nodal pts)
elx1 = (0:Nel1)*d1;
elx2 = (0:Nel2)*d2;

% Now calculate the collocation points using the Greville Abscissae

% To avoid repeats (since the first and last point are the same), only do
% n-1 pts for the first one. Following is # of col points
ncx1 = n1;
ncx2 = n2;

% Now for actual locations
cx1 = zeros(ncx1,1);
cx2 = zeros(ncx2,1);
for i = 1:ncx1
    cx1(i) = sum(t1(i+1:i+p1))/p1;
end
% cx1(end) = cx1(end)-.1*t1(end);

for i = 1:ncx2
    cx2(i) = sum(t2(i+1:i+p2))/p2;
end
% cx2(1) = 0.1*t2(end);
% cx2(end) = cx2(end)-0.1*t2(end);

% Now for the quadrature points.
ng1 = 40;
ng2 = 40;

% Locations:
gx1 = zeros(Nel1,ng1);
gx2 = zeros(Nel2,ng2);

% Corresponding locations in the parent domain for each element
for i = 1:Nel1
    [gx1(i,:),wg1] = lgwt(ng1,elx1(i+1),elx1(i));
end

for i = 1:Nel2
    [gx2(i,:),wg2] = lgwt(ng1,elx1(i+1),elx1(i));
end

wg1 = -wg1*2;
wg2 = -wg2*2;

% Now let's calculate the locations of the Gauss Points in physical space
% First get valuers of B-spline basis functions, N and M
[N, M] = BS(p1, p2, n1,n2, t1, t2, gx1, gx2);

% Now the NURBS basis functions
R = NURB(N,M,w);

% Now for actual physical coordinates
X = NURBev(R,P);

% Also need basis function for collocation points
[Nx, Mx] = BS(p1, p2, n1,n2, t1, t2, cx1, cx2);
Rx = NURB(Nx,Mx,w);

% Only a few shape functions (p+1 per element) are nonzero in a given
% element, so we won't need to loop over the ones that are zero later. They 
% will always be monotonically increasing, so we just need the first
% nonzero one for each element

% There's probably an analytical way to do this, but this is fast enough
% honestly
con1 = zeros(1,Nel1);
con2 = zeros(1,Nel2);
for i = 1:Nel1
    for j = 1:n1
        tmp = N{j,end};
        if tmp(i,1) ~= 0 
            con1(i) = j;
            break;
        end
    end
end

for i = 1:Nel2
    for j = 1:n2
        tmp = M{j,end};
        if tmp(i,1) ~= 0 
            con2(i) = j;
            break;
        end
    end
end

%% Pre-processing required each time step
% Calculate derivatives at each Gauss point
dN = BSd(N,t1,gx1);
dM = BSd(M,t2,gx2);
[dX1,dX2] = NURBd(N,M,dN,dM,P,w);

% Calculate Jacobian vector/magnitude and normal vector at Gauss points
% This is a cross product of dr/dx X dr/deta, where
% r is the parameterization of the splines.
[c,cm, n] = Ncross(dX1,dX2);

% Info to calculate integrand
xc = NURBev(Rx,P);

%% Calculate the integral
A = zeros(3*ncx1*ncx2,3*ncx1*ncx2);
b = zeros(3*ncx1*ncx2,1);

% Loop over collocation points
for ic1 = 1:ncx1
    for ic2 = 1:ncx2
%       What collocation point are we at?
        ind1 = ic2 + (ic1-1)*ncx2;
%       Column to put back into matrix
        row = 3*(ind1-1)+1;
        
%       Loop over elements
        for je1 = 1:Nel1
            for je2 = 1:Nel2
%               Loop over shape functions in element
                for kshp1 = con1(je1):con1(je1) + p1
                    for kshp2 = con2(je2):con1(je2) + p2
%                       Values of NURBS shape functions at current location
                        Rc = R{kshp1,kshp2};
%                       Loop over Gauss points at shape functions
                        At = zeros(3,3);
                        bt = zeros(3,1);
                        for lg1 = 1:ng1
                            for lg2 = 1:ng2
%                               Get x-y and normal vector for current combo
%                               of collocation pt and Gauss point, as well
%                               as Jacobian, Gauss weights, shape function
%                               value/location
                                r = xc(:,ic1,ic2) - X(:,je1,je2,lg1,lg2);
                                nc = n(:,je1,je2,lg1,lg2);
                                J = cm(je1,je2,lg1,lg2);
                                At = At + ...
                                    Rc(je1,je2,lg1,lg2)*Gij(r)*wg1(lg1)*wg2(lg2)*J*d1/2*d2/2;
%                               I'm lucky the below worked. In the future,
%                               we will want U to be the control points for
%                               the velocity, which we will need to
%                               calculate a priori. Only works cause U constant and partition of unity Or maybe we could just
%                               move this outside the shps loop?
                                bt = bt + ...
                                    Rc(je1,je2,lg1,lg2)*Tij(r,nc)*U*wg1(lg1)*wg2(lg2)*J*d1/2*d2/2;
                            end
                        end
                        
                        At = At/4/pi/mu;
                        bt = bt/4/pi;

%                       What shape functions are we integrating wrt to?
                        ind2 = kshp2 + (kshp1-1)*ncx2;
                        
%                       Get the column to put these into big mats
                        col = 3*(ind2-1)+1;
                        
%                       Put back together
                        b(row:row+2) = b(row:row+2) + bt;
                        
                        A(row:row+2,col:col+2) = ...
                            A(row:row+2,col:col+2) + At;
                    end %kshp2
                end %kshp1
            end %je2
        end %je1
        b(row:row+2) = b(row:row+2) - U;
    end
end

% Take care of repeated nodes
for ic1 = 1:ncx1
    for ic2 = 1:ncx2
        ind1 = ic2 + (ic1-1)*ncx2;
        row = 3*(ind1-1)+1;
        for khshp1 = 1:n1
            for kshp2 = 1:n2
                ind2 = kshp2 + (kshp1-1)*ncx2;
                col = 3*(ind2-1)+1;
                
%               We have wrapped all the way around set equal to beginning
                if (ic1 == ncx1)
                    A(row:row+2,:) = 0;
                    A(row,row) = 1;
                    A(row+1,row + 1) = 1;
                    A(row+2,row + 2) = 1;
                    A(row,row - 3*(ncx1-1)*(n2)) = -1;
                    A(row+1,row+1 - 3*(ncx1-1)*(n2)) = -1;
                    A(row+2,row+2 - 3*(ncx1-1)*n2) = -1;
                    b(row) = 0;
                    b(row+1) = 0;
                    b(row+2) = 0;
                end
                
                
%               Save column to set other endpoints to
                if (ic2 == n2 && ic1 == 1)
                    sacol1 = row;
                end
                
                if (ic2 == n2 && ic1 ~=1)
                    A(row:row+2,:) = 0;
                    A(row,row) = 1;
                    A(row+1,row + 1) = 1;
                    A(row+2,row + 2) = 1;
                    A(row,sacol1) = -1;
                    A(row+1,sacol1+1) = -1;
                    A(row+2,sacol1+2) = -1;
                    b(row) = 0;
                    b(row+1) = 0;
                    b(row+2) = 0;
                end
                
                
                if (ic2 == 1 && ic1 == 1)
                    sacol2 = row;
                end
                
                if (ic2 == 1 && ic1 ~=1)
                    A(row:row+2,:) = 0;
                    A(row,row) = 1;
                    A(row+1,row + 1) = 1;
                    A(row+2,row + 2) = 1;
                    A(row,sacol2) = -1;
                    A(row+1,sacol2+1) = -1;
                    A(row+2,sacol2+2) = -1;
                    b(row) = 0;
                    b(row+1) = 0;
                    b(row+2) = 0;
                end                
            end
        end
    end
end


F = A\b;
Ft = zeros(3,n1,n2);

% Putting it back together
cnt = 1;
for ic1 = 1:ncx1
    for ic2 = 1:ncx2
        Ft(1:3,ic1,ic2) = F(cnt:cnt + 2);
        cnt = cnt + 3;
    end
end
F = Ft;

% Reconstruct stresses at collocation points
% fc = NURBev(Rx,F);

%%
% % Plot tests

% % Some GOOD plots of pressures and whatnot at Gauss points
fc = NURBev(R,F);
cnt = 0;
cnt2 = 0;
clear x y z f1 f2 f3 xx ff nf
for i = 1:Nel1
    for j = 1:Nel2
        for k = 1:ng1
            for l =1:ng2
                cnt=cnt+1;
                nc(:,cnt) = n(:,i,j,k,l);
                if(i==3 && k==3)
                    cnt2 = cnt2+1;
                    xx(cnt2) = X(2,i,j,k,l);
                    ff(cnt2) = dot(fc(:,i,j,k,l),-nc(:,cnt));
                end
                x(cnt)=X(1,i,j,k,l);
                y(cnt)=X(2,i,j,k,l);
                z(cnt)=X(3,i,j,k,l);
                f1(cnt) = fc(1,i,j,k,l);
                f2(cnt) = fc(2,i,j,k,l);
                f3(cnt) = fc(3,i,j,k,l);
                nf(cnt) = dot(fc(:,i,j,k,l),-nc(:,cnt));%norm(fc(:,i,j,k,l));
            end
        end
    end
end
thet = acos(xx);
pa = -1.5*cos(thet);
% scatter3(x,y,z,100,nf,'filled')
hold on
% quiver3(x,y,z,f1,f2,f3)
% pbaspect([1,1,1])
% axis([-1,1,-1,1,-1,1])
% figure
plot(thet,ff,'o');
hold on
plot(thet,pa);
L2err = sqrt(sum(((ff-pa)./pa))^2)/length(ff);
avgerr = mean(abs(ff-pa)./abs(pa));
xlabel('\theta')
ylabel('P')
legend('BEM NURBS','Analytical','location','northwest')


% % Scatter of Gp
% x=[];
% y=[];
% z=[];
% for i = 1:Nel1
%     for j=1:Nel2
%         for k=1:ng1
%             for l= 1:ng2
%                 x(end+1)=X(1,i,j,k,l);
%                 y(end+1)=X(2,i,j,k,l);
%                 z(end+1)=X(3,i,j,k,l);
%             end
%         end
%     end
% end
% scatter3(x,y,z,100,'filled');
% pbaspect([1,1,1])
% hold on
% Scatter of collpts
% x2=[];
% y2=[];
% z2=[];
% for k=1:ncx1
%     for l= 1:ncx2
%         x2(end+1)=xc(1,k,l);
%         y2(end+1)=xc(2,k,l);
%         z2(end+1)=xc(3,k,l);
%     end
% end     
% scatter3(x2,y2,z2,100);
% pbaspect([1,1,1])
% hold on   

% % Scatter of Gp derivs
% for i = 1:Nel1
%     for j=1:Nel2
%         for k=1:ng1
% x=[];y=[];z=[];
% x2=[];y2=[];z2=[];
%             for l= 1:ng2
%                 x(end+1)=dX1(1,i,j,k,l);
%                 y(end+1)=dX1(2,i,j,k,l);
%                 z(end+1)=dX1(3,i,j,k,l);
%                 x2(end+1)=dX2(1,i,j,k,l);
%                 y2(end+1)=dX2(2,i,j,k,l);
%                 z2(end+1)=dX2(3,i,j,k,l);
%             end
% %         scatter3(x,y,z,100,'filled');
% %         hold on        
% %         scatter3(x2,y2,z2,100,'filled');
% %         pbaspect([1,1,1])
%         end
%     end
% end
% scatter3(x,y,z)

% 
% % Shape fns in first element
% 
% for i = 1:n1
%     tt = N{i,end};
%     tt = tt(1,:);
%     plot(gx1(1,:),tt)
%     hold on
% end
% figure
%
% % Derivs of shape fns in first element
% for i = 1:n1
%     ss = dN{i};
%     ss = ss(1,:);
%     plot(gx1(1,:),ss)
%     hold on
% end
% 
% % Shape fns in second element
% for i = 1:n1
%     tt = N{i,end};
%     tt = tt(2,:);
%     plot(gx1(2,:),tt)
%     hold on
% end
% % 
% % Normal vectors
% cnt = 0;
% nx = [];
% ny = [];
% nz = [];
% x = [];
% z=[];
% y = [];
% for i = 1:Nel1
%     for j = 1:Nel2
%         for k = 1:ng1
%             for l =1:ng2
%                 cnt=cnt+1;
%                 nx(cnt) = n(1,i,j,k,l);
%                 ny(cnt) = n(2,i,j,k,l);
%                 nz(cnt) = (3,i,j,k,l);
%                 x(cnt)=X(1,i,j,k,l);
%                 y(cnt)=X(2,i,j,k,l);
%                 z(cnt)=X(3,i,j,k,l);
%             end
%         end
%     end
% end
% fc = NURBev(Rx,F);
% cnt = 0;
% clear x y z f1 f2 f3
% for i = 1:ncx1
%     for j = 1:ncx2
%         cnt=cnt+1;
%         x(cnt)=xc(1,i,j);
%         y(cnt)=xc(2,i,j);
%         z(cnt)=xc(3,i,j);
%         f1(cnt) = fc(1,i,j);
%         f2(cnt) = fc(2,i,j);
%         f3(cnt) = fc(3,i,j);
%     end
% end
% 

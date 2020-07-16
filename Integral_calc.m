%% Paper in question
% https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
% stokes-flow-past-a-particle-of-arbitrary-shape-a-numerical-method-of-solution/
% 1DD984501FDD52C7F88811182CD4D286

%% Important!!
% Need to look into which collocation points (where to eval function) are
% most accurate. It seems that this is likely ON THE KNOT VECTOR knots??

%% Gauss points info
% Weights
gw = [2,0,0,0;
      1,1,0,0;
      5/9,8/9,5/9,0;
      (18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
    
% Locations
gloc = [0,0,0,0;
        -1/sqrt(3),1/sqrt(3),0,0;
        -sqrt(3/5),0,sqrt(3/5),0;
        -sqrt(3/7 + 2/7*sqrt(6/5)),-sqrt(3/7 - 2/7*sqrt(6/5)),sqrt(3/7 - 2/7*sqrt(6/5)),sqrt(3/7 + 2/7*sqrt(6/5))];

%% Actual code - NURBS Construction
% First lets get our stencil, which we need because the actual values of xi
% might vary in the knot vector, especially with periodic NURBS
% We will only use certain knot spans (Again depending on if this is
% periodic or not) and the stencil just tells us the points at various
% proportions along that knot span we will use. Basically a normalization
% of the xi coordinates

% Number of surfaces in each direction as well as the width

% Very important!!! This is the most accurate when the surfaces are
% approximately square. I don't know why, but Sxi = Seta = 100 is much less
% accurate than Sxi = 40 and Seta = 20
Sxi = 8;
Seta = 4;

dxi = 1/Sxi;
deta = 1/Seta;

% Number of Gauss points on a surface in each direction and Gauss weights
ngxi = 4;
ngeta = 4;

% Now we know how many points in our stencil, e.g. each surface needs ngxi
% x ngeta points to numerically integrate
xis = zeros(ngxi*Sxi,1);
etas = zeros(ngeta*Seta,1);
wgxi = zeros(1,ngxi);
wgeta= zeros(1,ngeta);

% Let's allocate the weights and locations of the xis/etas in the first
% surface
for i = 1:ngxi
    wgxi(i) = gw(ngxi,i);
    xis(i) = dxi/2 + (dxi/2*gloc(ngxi,i));
end

for i = 1:ngeta
    wgeta(i) = gw(ngeta,i);
    etas(i) = deta/2 + (deta/2*gloc(ngeta,i));
end

% Now for the rest of the points just add d -xi/-eta to the first points
for i = 2:Sxi
    for j = 1:ngxi
        xis((i-1)*ngxi + j) = xis((i-2)*ngxi + j) + dxi;
    end
end

for i = 2:Seta
    for j = 1:ngeta
        etas((i-1)*ngeta + j) = etas((i-2)*ngeta + j) + deta;
    end
end

% Actual NURBS construction now that we know points.

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
t1 = [0,0,0,.25,.25,.5,.5,.75,.75,1,1,1];
t2 = [0,0,0,0.5,0.5,1,1,1];

% Control net

P = zeros(4,n1,n2);
xs = [1,1,0,-1,-1];
ys  = [0,1,1,1,0];
r2 = sqrt(2)/2;
w1 = [1,r2,1,r2,1,r2,1,r2,1];
w2 = [1,r2,1,r2,1];

for i = 1:n1
    P(2,i,:) = xs;
    P(4,i,:) = w2*w1(i);
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

% Evaluate basis functions and points at stencils
[N, M] = BS(p1, p2, n1,n2, t1, t2, xis, etas);
X = NURB(N,M,P);

% Now derivatives at stencils and points
dM = BSd(M,t2,etas);
dN = BSd(N,t1,xis);
[Xx,Xe] = NURBd(N,M,dN,dM,P);

% Cross product results. This is a cross product of dr/dx X dr/deta, where
% r is the parameterization of the splines.
[c,cm, n] = Ncross(Xx,Xe);

%% Info to calculate integrand
% Need xm, the centers of the surfaces where we're evaluating these things
% Just take these as halfway through the surfaces in parent domain
% Follows similar process as above
xmxis  = linspace(dxi /2,1-dxi /2,Sxi);
xmetas = linspace(deta/2,1-deta/2,Seta);
[Nx, Mx] = BS(p1, p2, n1,n2, t1, t2, xmxis, xmetas);
xm = NURB(Nx,Mx,P);

U = [0,1,0];

%% Calculate the integral
% To save time, I think that eventually I will need to do the big integral
% that I just break up into the sums as an acutal big integral
g = zeros(3*Sxi*Seta,1);
A = zeros(3*Sxi*Seta,3*Sxi*Seta);

%Loop over all the surfaces in each direction
% Loop over xm surfaces
for i =1:Sxi
    for j = 1:Seta
%       Dimension loop
        for d1 = 1:3
%           Current index
            ii = d1 + (j-1)*3 + (i-1)*Seta*3;
%           Loop over y surfaces
%           Total integral for the given xm
            I = 0;
            for i2 = 1:Sxi
                for j2 = 1:Seta
%                   A integrals
                    for d2 = 1:3
%                       Current index
                        jj = d2 + (j2-1)*3 + (i2-1)*Seta*3;
%                       Loop over Gauss points
                        SI = 0;
                        for k = 1:ngxi
                            for l = 1:ngeta
                                ixi = (i2-1)*ngxi + k;
                                jeta= (j2-1)*ngeta+ l;
                                r = sqrt((xm(1,i,j)-X(1,ixi,jeta))^2 ...
                                        +(xm(2,i,j)-X(2,ixi,jeta))^2 ...
                                        +(xm(3,i,j)-X(3,ixi,jeta))^2);
                                Atmp = (xm(d1,i,j) - X(d1,ixi,jeta)) ...
                                     * (xm(d2,i,j) - X(d2,ixi,jeta)) ...
                                     / r^3;
                                 if d1 == d2
                                     Atmp = Atmp +1/r;
                                 end
                                 SI = SI + Atmp*cm(ixi,jeta)*wgxi(k)*wgeta(l)*(deta/2)*(dxi/2);
                            end
                        end
                        A(ii,jj) = SI/4/pi;
                    end
                    
%                   g integrals
%                   Gauss points in y surfaces
%                   Just the integral for the one surface
                    SI = 0;
                    for k = 1:ngxi
                        for l = 1:ngeta
                            ixi = (i2-1)*ngxi + k;
                            jeta= (j2-1)*ngeta+ l;
                            r = sqrt((xm(1,i,j)-X(1,ixi,jeta))^2 ...
                                    +(xm(2,i,j)-X(2,ixi,jeta))^2 ...
                                    +(xm(3,i,j)-X(3,ixi,jeta))^2);
                             
                            gtmp =((xm(d1,i,j) - X(d1,ixi,jeta)) ...
                                 * ( (xm(1,i,j) - X(1,ixi,jeta))*n(1,ixi,jeta) ...
                                 +   (xm(2,i,j) - X(2,ixi,jeta))*n(2,ixi,jeta) ...
                                 +   (xm(3,i,j) - X(3,ixi,jeta))*n(3,ixi,jeta))...
                                 * ( (xm(1,i,j) - X(1,ixi,jeta))*U(1) ...
                                 +   (xm(2,i,j) - X(2,ixi,jeta))*U(2) ...
                                 +   (xm(3,i,j) - X(3,ixi,jeta))*U(3)))/r^5;
                             SI = SI + gtmp*cm(ixi,jeta)*wgxi(k)*wgeta(l)*(deta/2)*(dxi/2);
                        end
                    end
                    I = I + SI;
                end
            end
            g(ii) = -U(d1) - 3/2/pi*I;
        end
    end
end
% A=A/2;
% Now we can calculate f
f = A\g;

%% Field reconstruction
% % Loop over xm surfaces
% 
% nx = 30;
% ny = 30;
% XX = linspace(-2,2,nx);
% YY = linspace(-2,2,ny);
% u = zeros(3,nx,ny);
% SI = zeros(1,3);
% I = SI;
% 
% for xs = 1:nx
% for ys = 1:ny
% % point to evaluate at
% myx = [XX(xs),YY(ys),0];
% I(:) = 0;
% for i =1:Sxi
%     for j = 1:Seta
%         ii = (j-1)*3 + (i-1)*Seta*3;
%         for d1 = 1:3
%             SI(:) = 0;
%             SI2 = 0;
%             for k = 1:ngxi
%                 for l = 1:ngeta
%                     ixi = (i-1)*ngxi + k;
%                     jeta= (j-1)*ngeta+ l;
%                     r = sqrt((myx(1)-X(1,ixi,jeta))^2 ...
%                             +(myx(2)-X(2,ixi,jeta))^2 ...
%                             +(myx(3)-X(3,ixi,jeta))^2);
%                         
% %                   First integral
%                     for d2 = 1:3
%                         utmp =((myx(d1) - X(d1,ixi,jeta)) ... 
%                              * (myx(d2) - X(d2,ixi,jeta)))/r^3;
%                         if d1 == d2
%                             utmp = utmp +1/r;
%                         end
%                         utmp = utmp/(-8*pi); 
%                         SI(d2) = SI(d2) + utmp*cm(ixi,jeta)*wgxi(k)*wgeta(l)*(deta/2)*(dxi/2);
%                     end
% 
% %                       Second Integral
%                     utmp = ((myx(d1) - X(d1,ixi,jeta)) ...
%                              * ( (myx(1) - X(1,ixi,jeta))*n(1,ixi,jeta) ...
%                              +   (myx(2) - X(2,ixi,jeta))*n(2,ixi,jeta) ...
%                              +   (myx(3) - X(3,ixi,jeta))*n(3,ixi,jeta))...
%                              * ( (myx(1) - X(1,ixi,jeta))*U(1) ...
%                              +   (myx(2) - X(2,ixi,jeta))*U(2) ...
%                              +   (myx(3) - X(3,ixi,jeta))*U(3)))/r^5;
%                     SI2 = SI2 + utmp*cm(ixi,jeta)*wgxi(k)*wgeta(l)*(deta/2)*(dxi/2);
%                 end
%             end
%             I(d1) = I(d1) ...
%                   - (SI(1)*f(ii+1) + SI(2)*f(ii+2) + SI(3)*f(ii+3)) ...
%                   + 3/4/pi*SI2;
%         end  
%     end
% end
% if (myx(1)^2 + myx(2)^2 + myx(3)^2) > 1.0
%     u(:,xs,ys) = I;
% else
%     u(:,xs,ys) = -U;
% end
% 
% end
% end
% 
% 
% UU = u;
% UU(1,:,:) = UU(1,:,:)+1;
% utt = zeros(nx,ny);
% vtt = zeros(nx,ny);
% for i = 1:nx
%     for j = 1:ny
%         utt(i,j)=u(1,i,j)+1;
%         vtt(i,j)=u(2,i,j);
%     end
% end

%% plots
% starty = linspace(min(YY), max(YY),nx);
% startx = ones(size(starty))*min(YY);
% quiver(XX,YY,utt,vtt)
% theta = linspace(0,2*pi,100);
% x = cos(theta);
% y = sin(theta);
% hold on
% plot(x,y);
% pbaspect([1,1,1])
% axis([-2,2,-2,2])

%% Let's just get the drag. This should be easy:
% f is just a surface stress vector, so we just need to sum up f_i*A_i, the
% the stress on an element times its area. I think the hardest part will be
% calculating the area
F = zeros(3,1);
Atmp = 0;

% loop over all the areas
for i =1:Sxi
    for j = 1:Seta
        ii = (j-1)*3 + (i-1)*Seta*3;
%       get area
%       But first get index to average cm through
        ixi = (i-1)*ngxi + 1;
        jeta = (j-1)*ngeta + 1;
        Atmp = mean(mean(cm(ixi+ngxi-1,jeta+ngeta-1)))*deta*dxi;
        F(1) = F(1) + Atmp*f(ii+1);
        F(2) = F(2) + Atmp*f(ii+2);
        F(3) = F(3) + Atmp*f(ii+3);
        
    end
end

%% Getting pressures
% Should just be normal portion of f
p = zeros(Sxi,Seta);
tau = p;

% Now derivatives at stencils and points
dMx = BSd(Mx,t2,xmetas);
dNx = BSd(Nx,t1,xmxis);
[Xxx,Xxe] = NURBd(Nx,Mx,dNx,dMx,P);

% Cross product results
[cx,cmx, nx] = Ncross(Xxx,Xxe);

for i = 1:Sxi
    for j = 1:Seta
        for d = 1:3
            ii = d + (j-1)*3 + (i-1)*Seta*3;
            p(i,j) = p(i,j) + f(ii)*nx(d,i,j);
        end
        tau(i,j) = sqrt(sqrt(f(ii-2)^2 + f(ii-1)^2 + f(ii)^2)^2-p(i,j)^2);
    end
end

% get points to plot along
xpl = zeros(1,Seta);
for i = 1:Seta
    xpl(i) = xm(2,1,i);
end
thet = acos(xpl);
ppl = p(1,:);
pa = 1.5*cos(thet);
taup= tau(1,:);
taua=1.5*sin(thet);
% figure
plot(thet,ppl);
hold on
plot(thet,pa);
plot(thet,taup);
plot(thet,taua);
hold off

% 1st order as of right now
L2err = sqrt(sum(((ppl-pa)./pa))^2)/length(ppl);
L1err = max(abs((ppl-pa)./pa));
avgerr= mean(abs((ppl-pa)./pa));

% axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
% pbaspect([1,1,1]);


%% numerical integration tests
% clear fun mid;
% nn=10;
% x = [0.5,0.5];
% y = linspace(0,.5,nn);
% [X,Y] = meshgrid(y,y);
% for i =1:nn
%     for j=1:nn
%         d(1) = x(1)-y(i);
%         d(2) = x(2)-y(j);
%         r = sqrt(d(1)^2 + d(2)^2);
%         fun(i,j) = d(1)*d(2)/r^3;
%         if isnan(fun(i,j))
%             fun(i,j) = fun(i-1,j);
%         end
%         if i>1&&j>1
%             mid(i,j) = 1/4*(fun(i,j)+fun(i-1,j)+fun(i,j-1)+fun(i-1,j-1));
%         end
%     end
% end
% dA=1/nn/nn;
% % intg=sum(sum(mid))*dA;
% % ns = nn;
% intg(end+1)=sum(sum(mid))*dA;
% ns(end+1) = nn;
% % disp(intg);

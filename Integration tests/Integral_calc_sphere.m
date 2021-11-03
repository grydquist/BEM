%% Gauss points info
% weights
gw = [2,0,0,0;
      1,1,0,0;
      5/9,8/9,5/9,0;
      (18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
    
% Locations
gloc = [0,0,0,0;
        -1/sqrt(3),1/sqrt(3),0,0;
        -sqrt(3/5),0,sqrt(3/5),0;
        -sqrt(3/7 + 2/7*sqrt(6/5)),-sqrt(3/7 - 2/7*sqrt(6/5)),sqrt(3/7 - 2/7*sqrt(6/5)),sqrt(3/7 + 2/7*sqrt(6/5))];

%% Actual code
% First lets get our stencil, which we need because the actual values of xi
% might vary in the knot vector, especially with periodic NURBS

% Number of surfaces in each direction as well as change in variables
Sxi =20;
Seta = 10;

dxi = 1/Sxi;
deta = 1/Seta;

% number of Gauss points on a surface in each direction and Gauss weights
ngxi = 3;
ngeta = 3;

% now we know how many points in our stencil
xis = zeros(ngxi*Sxi,1);
etas = zeros(ngeta*Seta,1);
wgxi = zeros(1,ngxi);
wgeta= zeros(1,ngeta);

% Let's allocate the weights and first locations of the xis/etas we want
for i = 1:ngxi
    wgxi(i) = gw(ngxi,i);
    xis(i) = dxi/2 + (dxi/2*gloc(ngxi,i));
end

for i = 1:ngeta
    wgeta(i) = gw(ngeta,i);
    etas(i) = deta/2 + (deta/2*gloc(ngeta,i));
end

% Now for the rest of the points just add dxi/eta to the first points
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

% Cross product results
[c,cm,n] = Ncross(Xx,Xe);

%% Calculate the integral
%Loop over all the surfaces in each direction
I = 0;
for i =1:Sxi
    for j = 1:Seta
        
        SI = 0;
%       Loop over Gauss points in this surface
        for k = 1:ngxi
            for l = 1:ngeta
                ixi = ngxi* (i-1) + k;
                ieta= ngeta*(j-1) + l;
                SI = SI + cm(ixi,ieta)*wgxi(k)*wgeta(l);
            end
        end
        SI = SI*(deta/2)*(dxi/2);
        I = I + SI;
    end
end

disp(I - 4*pi);

% figure
% scatter3(X(1,:),X(2,:),X(3,:));
% axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
% pbaspect([1,1,1]);
% 
% figure
% for i = 1:n1
%     plot(xi,N{i,end})
%     hold on
% end
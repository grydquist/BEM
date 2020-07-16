%% Evaluate the 1st set of basis functions
% order
p1 = 3;
% number of control points/basis functions
n1 = 7;
% size of knot vector
k1 = n1+p1+1;
% actual knot vector
t1 = linspace(0,1,k1);

% iter is number of basis functions to evaluate for each N
iter = n1+p1;
% N contains all basis functions of all orders
N = cell(n1 + p1,p1);
% Evaluation points
ev1 = 100;
% tmp is the temporary holder for the evaluation of a given basis function
tmp = zeros(1,ev1);
% basis functions of 1 order lower at i and i+1
tmpm=tmp;
tmpp=tmp;
% xi is the knot value
xi = linspace(min(t1),max(t1),ev1);

% Outer loop goes from order 0 -> order p basis functions
for i = 1:p1+1
%   Inner loop evaluates the basis functions
    for j = 1:iter
%       Get lower order basis functions
        if i>1
            tmpm = N{j,i-1};
            tmpp = N{j+1,i-1};
        end
%       Loop over knot values to eval
        for kk = 1:ev1
%           Evaluate 0th order basis fns
            if i==1
                if xi(kk)< t1(j+1)&& xi(kk)>=t1(j)
                    tmp(kk) = 1;
                else
                    tmp(kk) = 0;
                end
%           Evaluate the ith order basis function    
            else
                tmp(kk) = 0;
%               This is to deal with repeated knots
                if t1(j+i-1)~=t1(j)
                    tmp(kk) = (xi(kk) - t1(j)  )/(t1(j+i-1) - t1(j)  )*tmpm(kk);
                end
                if t1(j+i) ~= t1(j+1)
                    tmp(kk) = tmp(kk) + (t1(j+i) - xi(kk))/(t1(j+i  ) - t1(j+1))*tmpp(kk);
                end
            end
        end
        N{j,i} = tmp;
    end
    iter = iter-1;
end

%% Now for the 2nd set of basis fns
% order
p2 = 2;
% number of control points/basis functions
n2 = 8;
% size of knot vector
k2 = n2+p2+1;
% actual knot vector
t2 = [0,0,0,1,2,3,4,5,6,6,6];

% iter is number of basis functions to evaluate for each M
iter = n2+p2;
% N contains all basis functions of all orders
M = cell(n2 + p2,p2);
% tmp is the temporary holder for the evaluation of a given basis function
ev2 = 50;
tmp = zeros(1,ev2);
% basis functions of 1 order lower at i and i+1
tmpm=tmp;
tmpp=tmp;
% nu is the knot value
nu = linspace(min(t2),max(t2),ev2);


% Outer loop goes from order 0 -> order P basis functions
for i = 1:p2+1
%   Inner loop evaluates the basis functions
    for j = 1:iter
%       Get lower order basis functions
        if i>1
            tmpm = M{j,i-1};
            tmpp = M{j+1,i-1};
        end
%       Loop over knot values to eval
        for kk = 1:ev2
%           Evaluate 0th order basis fns
            if i==1
                if nu(kk)< t2(j+1)&& nu(kk)>=t2(j)
                    tmp(kk) = 1;
                else
                    tmp(kk) = 0;
                end
%           Evaluate the ith order basis function    
            else
                tmp(kk) = 0;
%               This is to deal with repeated knots
                if t2(j+i-1)~=t2(j)
                    tmp(kk) = (nu(kk) - t2(j)  )/(t2(j+i-1) - t2(j)  )*tmpm(kk);
                end
                if t2(j+i) ~= t2(j+1)
                    tmp(kk) = tmp(kk) + (t2(j+i) - nu(kk))/(t2(j+i  ) - t2(j+1))*tmpp(kk);
                end
            end
        end
        M{j,i} = tmp;
    end
    iter = iter-1;
end

%% Control net and weighted denominator
P = zeros(4,n1,n2);
xzs = [0,.3,.58,1,1,.58,.3,0];
ys  = [.1033,.1033,.3270,.3270,-.3270,-.3270,-.1033,-.1033];
w =[1,0.3,1,0.3,0.3,1,0.3,1];
P(1,1,:) = zeros(1,n2)';
P(1,2,:) = -1.5*xzs;%-2*xzs;
P(1,3,:) = 0;
P(1,4,:) = 1.5*xzs;%1*xzs;
P(1,5:7,:) = P(1,1:3,:);

for i = 1:n1
    P(2,i,:) = ys;
    P(4,i,:) = w;
end

P(3,1,:) = 1.5*xzs;%0.75*xzs;
P(3,2,:) = 0;
P(3,3,:) = -1.5*xzs;%-2.0*xzs; %%%% pulling here 1.5 orig
P(3,4,:) = 0;
P(3,5:7,:) = P(3,1:3,:);



NMw = zeros(ev1,ev2);
for i = 1:n1
    tmp1 = N{i,end};
    for j = 1:n2
        tmp2 = M{j,end};
        for k = 1:ev1
            for l = 1:ev2
                NMw(k,l) = NMw(k,l) + tmp1(k)*tmp2(l)*P(4,i,j);
            end
        end
    end
end

%% Evaluate at a bunch of points
uu = zeros(ev1,ev2);
vv = uu;
ww = vv;
track = [];

% Loop over first set of basis functions
for i =1:n1
%   ith basis function of N
    tmp1 = N{i,end};
    for j = 1:n2
%       jth basis function of M
        tmp2 = M{j,end};
%       Loop over evaluations on nu and xi
        for k = 1:ev1
            for l = 1:ev2
%               Evaluate only if in these knot spans
                if xi(k)>=t1(p1+1)&&xi(k)<=t1(n1+1) && nu(l)>=t2(p2+1)&&nu(l)<=t2(n2+1)
                    uu(k,l) = uu(k,l) + tmp1(k)*tmp2(l)*P(1,i,j)*P(4,i,j)/NMw(k,l);
                    vv(k,l) = vv(k,l) + tmp1(k)*tmp2(l)*P(2,i,j)*P(4,i,j)/NMw(k,l);
                    ww(k,l) = ww(k,l) + tmp1(k)*tmp2(l)*P(3,i,j)*P(4,i,j)/NMw(k,l);
                end
            end            
        end
    end
end

%% Plot points
% reshape into vectors
up  =reshape(uu,ev1*ev2,1);
vp  =reshape(vv,ev1*ev2,1);
wp  =reshape(ww,ev1*ev2,1);

% remove points
up(track) = [];
vp(track) = [];
wp(track) = [];
% figure
% scatter3(up,vp,wp);
% axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
% pbaspect([1,1,1]);
%% plot basis functions
% for i = 1:n2
%     plot(eta,M{i,end})
%     hold on
% end
% axis([-.1,1.1,-.1,1.1])
plot(up)

% order
p = 3;
% number of control points/basis functions
n = 7;
% size of knot vector
k = n+p+1;
% actual knot vector
t = linspace(0,1,k);
% actual control points
P = zeros(2,n);

% circle
a = 2*pi/(n-p);
h = sin(a)/(2+cos(a));
mu = 3/(2+cos(a));
 
for i = 1:(n-p)
   P(1,i) = cos((i)*a)*mu;
   P(2,i) = sin((i)*a)*mu;
end
w = ones(1,n);

P(1:2,end-(p-1):end) = P(1:2,1:p);
w(1,end-(p-1):end) = w(1,1:p);

% % Nurbs circle
% P = [1,1,0,-1,-1,-1,0,1,1;
%      0,1,1,1,0,-1,-1,-1,0];
% t = [0,0,0,1,1,2,2,3,3,4,4,4];
% ww = sqrt(2)/2;
% w = [1,ww,1,ww,1,ww,1,ww,1];

%% Now let's evaluate the basis functions
% iter is number of basis functions to evaluate for each N
iter = n+p;
% N contains all basis functions of all orders
N = cell(n + p,p);
% Evaluation points
ev = 1000;
% tmp is the temporary holder for the evaluation of a given basis function
tmp = zeros(1,ev);
% basis functions of 1 order lower at i and i+1
tmpm=tmp;
tmpp=tmp;
% xi is the knot value
xi = linspace(min(t),max(t),ev);

% Outer loop goes from order 0 -> order P basis functions
for i = 1:p+1
%   Inner loop evaluates the basis functions
    for j = 1:iter
%       Get lower order basis functions
        if i>1
            tmpm = N{j,i-1};
            tmpp = N{j+1,i-1};
        end
%       Loop over knot values to eval
        for kk = 1:ev
%           Evaluate 0th order basis fns
            if i==1
                if xi(kk)<= t(j+1)&& xi(kk)>=t(j)
                    tmp(kk) = 1;
                else
                    tmp(kk) = 0;
                end
%           Evaluate the ith order basis function    
            else
                tmp(kk) = 0;
%               This is to deal with repeated knots
                if t(j+i-1)~=t(j)
                    tmp(kk) = (xi(kk) - t(j)  )/(t(j+i-1) - t(j)  )*tmpm(kk);
                end
                if t(j+i) ~= t(j+1)
                    tmp(kk) = tmp(kk) + (t(j+i) - xi(kk))/(t(j+i  ) - t(j+1))*tmpp(kk);
                end
            end
        end
        N{j,i} = tmp;
    end
    iter = iter-1;
end
% Get the weighted basis fns
R = cell(1,n);
Nw = zeros(1,ev);

% Denominator (2.27 in ch. 2)
for i = 1:n
    Nw = Nw + N{i,end}*w(i);
end

for i = 1:n
    R{i} = N{i,end}*w(i)./Nw;
end

% We only want to evaluate these functions in the knot spans between p+1
% and n+1

% Evaluate using control points
uu = zeros(1,ev);
vv = uu;
track = [];
for i = 1:n
    tmp = R{i};
    for j = 1:ev
%       Only eval on above knot spans
        if xi(j)>t(p+1)&&xi(j)<t(n+1)
            uu(j) = uu(j) + tmp(j)*P(1,i);
            vv(j) = vv(j) + tmp(j)*P(2,i);
%       Keep track of these spots for removal
        else
            track = horzcat(track,j);
        end
    end
end
% Remove unused points
uu(track) = [];
vv(track) = [];
% Close loop
uu = horzcat(uu,uu(1));
vv = horzcat(vv,vv(1));

hold on;
plot(uu,vv)
axis([-1.5,1.5,-1.5,1.5])
pbaspect([1,1,1]);
%% plot basis functions
% % for i = 1:n
% %     plot(xi,N{i,end})
% %     hold on
% % end
% % axis([-.1,1.1,-.1,1.1])

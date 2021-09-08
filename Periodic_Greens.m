%% Calculates the periodic green's function for just one point,
% arbitrary lattice
clear eps1
global eps1

% First, define the lattice (columns vectors)
bv = [1 0 0; 0 1 0; 0 0 1];
x0 = [.5;.5;.5];
x = x0;
xh0 = x-x0;

% reciprocal basis vectors
kv = bv;
tau = dot(cross(bv(:,1),bv(:,2)),bv(:,3));
kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3));
kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1));
kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2));

eps1 = 0001;% pi^0.5/tau^(1/3);

%% Calculate them!
G = zeros(3,3);
G2 = G;
tst = zeros(41*41*41-1,1);
tst2 = tst;
r = tst;
ks = r;
% loop over 20 boxes in each direction
it = 0;
% for i = -20:20
%     for j = -20:20
%         for k = -20:20
%             if (i==0 && j==0 && k==0);continue;end
%             it = it + 1;
%             rc = i*bv(:,1) + j*bv(:,2) + k*bv(:,3);
%             xn = x + rc - x0;
%             kc = i*kv(:,1) + j*kv(:,2) + k*kv(:,3);
%             G = G + shortG(xn);
%             if ~(i==0 && j==0 && k==0)
%                 G = G + specG(kc)/tau*exp(-norm(kc)^2/4*eps1^2);
% %                 tst(it) = norm(specG(kc,xh0))/tau;
%             end
%             tmp = shortG(xn) + specRealG(xn) - Gij(xn);
%             tst(it) = norm(tmp);
% %             tst(it) = norm(shortG(xn));
%             r(it) = norm(xn);
%             ks(it) = norm(kc);
%             G2 = G2 + Gij(xn);
%             tmp = Gij(xn);
%             tst2(it) = tmp(1);
%             if(r(it)<1)
% %                 disp('test')
%             end
%         end
%     end
% end

%% Try ust calculating velocity first?? Try
u = [0;0;0];
u2 = u;
f = [1;0;0];
x0 = [.5;.5;.5];
V = 1;
it = 0;
uSR = zeros(41*41*41-1,1);
uSM = uSR;
ur = uSM;
r = uSR;
for i = -20:20
    for j = -20:20
        for k = -20:20
            if (i==0 && j==0 && k==0);continue;end
            it = it+1;
%           Vectors (real and reciprocal)            
            p = [i;j;k];
            kv= 2*pi*[i;j;k];
            
%           Real part
            A = shortG(p);
            if A(1,1)<-.05
                disp('?')
            end
%           Reciprocal part
            B = specG(kv);
            
%           Add in contributions
            tmp =  A*f;
            uSR(it) = tmp(1);
            u = u + A*f;
             
            tmp = exp(-norm(kv)^2/4/eps1^2)/V*B*f*exp(-1i*dot(kv,x0-x0));
            uSM(it) = tmp(1);
            u = u + exp(-norm(kv)^2/4/eps1^2)/V*B*f*exp(-1i*dot(kv,x0-x0));
            
            r(it) = norm(p);
%           Normal sum??
            tmp = Gij(p)*f;
            ur(it) = tmp(1);
            u2 = u2 + Gij(p)*f;
            
        end
    end
end
% Self-interaction
us = 8*eps1/sqrt(pi)*f;
% u = u -us;

%% Do for the Green's function 2/r
G = 0;
G2= 0;
x0 = [.5;.5;.5];
V = 1;
it = 0;
GSR = zeros(41*41*41-1,1);
GSM = GSR;
Gr = GSM;
r = GSR;
eps1 = 0.1;
for i = -20:20
    for j = -20:20
        for k = -20:20
            if (i==0 && j==0 && k==0);continue;end
            it = it+1;
            p = [i;j;k];
            kv= [i;j;k]*2*pi;
            
            G = G + rspecG(kv)/V;
            G = G + rshortG(p);
            G2 = G2 + 2/norm(p);
            r(it) = norm(p);
            GSM(it) = rspecG(kv)/V;
            GSR(it) = rshortG(p);
            Gr(it) = 2/norm(p);
        end
    end
end

%% Green's function functions
function G = rshortG(x)
global eps1
r = norm(x);
er = eps1*r;
G = 2/r*(2*er/sqrt(pi)*exp(-er^2)*(er^2-2) + erfc(er));
end

function G = rspecG(k)
global eps1
kn = norm(k);
w = kn/eps1;
G = 8*pi/kn^2*(1/w^4 + 1/4/w^2 + 1/8)*exp(-w^2/4)*1; % 1 for now since no sum
end

function G = shortG(x)
global eps1
r = norm(x);
er = r*eps1;
fp = erfc(er) - 2*er*exp(-er^2)/sqrt(pi);
fppr= -4*er*exp(-er^2)/sqrt(pi) + 4*er^3*exp(-er^2)/sqrt(pi);
C = fppr + fp;
D = fp - fppr;
G = eye(3)*C/r + x*x'*D/r^3;
end

% function G = shortG(x)
% global eps1
% r = norm(x);
% er = r*eps1;
% C = erfc(er) - 2/sqrt(pi)*(2*er^2 - 3)*er*exp(-er^2);
% D = erfc(er) + 2/sqrt(pi)*(1 - 2*er^2)*er*exp(-er^2);
% G = eye(3)*C/r + x*x'*D/r^3;
% 
% end

function G = specRealG(x)
global eps1
r = norm(x);
er = r*eps1;
fp = erf(er) + 2*er*exp(-er^2)/sqrt(pi);
fppr= 4*er*exp(-er^2)/sqrt(pi) - 4*er^3*exp(-er^2)/sqrt(pi);
C = fppr + fp;
D = fp - fppr;
G = eye(3)*C/r + x*x'*D/r^3;
end

function G = specG(k)
global eps1
kn = norm(k);
w = kn/eps1;
G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2 + 1/8)*(kn^2*eye(3) - k*k'); % No 1/k^2???
end

% function Gh = specG(k,x0) %%% STRANGE STUFF HERE
% global eps1
% w = norm(k)/eps1;
% kn = norm(k);
% Gh = 8*pi/(kn^2)*(eye(3) - k*k'/kn^2)*(1 + 1/4*w^2 + 1/8*w^4)*exp(-1/4*w^2)*exp(1i*dot(k,x0));
% end

%% Hasimotos
function G = HspecG(k)
global eps1
kn = norm(k);
G = 8*pi*(1 + kn^2/4/eps1^2)/kn^4*(kn^2*eye(3) - k*k');
end

function G = HshortG(x)
global eps1
r = norm(x);
er = eps1*r;
G = 2*(eps1*exp(-er^2)/sqrt(pi)/r^2 + erfc(er)/2/r^2)*(r^2*eye(3) + x*x') - 4*eps1/sqrt(pi)*exp(-er^2)*eye(3);
end


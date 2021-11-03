%% Calculates the periodic green's function for just one point,
% arbitrary lattice - not yet
clear eps1
global eps1

% First, define the lattice (columns vectors)
bv = [1 0 0; 0 1 0; 0 0 1];
x0 = [.25,.5,.5
      .75,.5,.5]';
x = [.5;.75,;.5];
xh0 = x0;
xh0(:,1) = x-x0(:,1);
xh0(:,2) = x-x0(:,2);

% reciprocal basis vectors
kv = bv;
tau = dot(cross(bv(:,1),bv(:,2)),bv(:,3));
kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3));
kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1));
kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2));

eps1 = 100;%pi^0.5/tau^(1/3); %00010.000;% 

f= [1,0,0
   -1,0,0]';

%% Calculate them! Calculate force at one point between two forces
it = 0;
bxs = 10;
bxt = bxs;
ucv = zeros(3,bxt);
uu=[0;0;0];
Gt = zeros(bxt,bxt,bxt);
Gt2 = Gt;
% Loop to check convergence with number  of boxes

for bxs = bxt:bxt
G = zeros(3,3,2);
Gh = G;
Gr = G;
% Loop over point forces
for ps = 1:2
for i = -bxs:bxs
    for j = -bxs:bxs
        for k = -bxs:bxs
            it = it + 1;
            
%           First we calculate the real space part
%           Get current distance and use with real space 
            rc = i*bv(:,1) + j*bv(:,2) + k*bv(:,3);
            xn = (x - x0(:,ps)) + rc;
            
%           Doing both Poz and Hasimoto
            G(:,:,ps)  = G(:,:,ps)  +  shortG(xn);
            Gh(:,:,ps) = Gh(:,:,ps) + HshortG(xn);
            
%           Now get the Fourier wave number associated with this box
%           Use that to get the Fourier part
            kc = i*kv(:,1) + j*kv(:,2) + k*kv(:,3);
            if(~(i==0 && j==0 && k==0))
                G(:,:,ps)  = G(:,:,ps)  +  specG(kc)/tau*exp(-1i*dot(kc,xh0(:,ps))); % do cosines??
                Gh(:,:,ps) = Gh(:,:,ps) + HspecG(kc)/tau*exp(-1i*dot(kc,xh0(:,ps)));
                
%               Just doing 3d fourier transform on real spectral
%                 tmp = specG(kc)/tau*exp(-1i*dot(kc,xh0(:,ps)));
%                 Gt2(i+bxt+1,j+bxt+1,k+bxt+1)= tmp(1,1);
            end
            
%           Also doing the real space one for comparison
            Gr(:,:,ps) = Gr(:,:,ps) + eye(3)/norm(xn) + xn*xn'/norm(xn)^3;
            
%           To check consistency with real space one            
%             G(:,:,ps)  = G(:,:,ps)  +  shortG(xn);
%             G(:,:,ps)  = G(:,:,ps)  +  specRealG(xn);

%           Just doing 3d fourier transform on real spectral
%             tmp = specG(kc)/tau*exp(-1i*dot(kc,xh0(:,ps)));
%             Gt(i+bxt+1,j+bxt+1,k+bxt+1) = tmp(1,1);

%           Just calc vel (real than fourier)
            uu = uu
            
        end
    end
end
end
% Only if particle is in the middle...
u = G(:,:,1)*f(:,1) + G(:,:,2)*f(:,2);
ucv(1,bxs) = real(u(2));

u = Gh(:,:,1)*f(:,1) + Gh(:,:,2)*f(:,2);
ucv(2,bxs) = real(u(2));

u = Gr(:,:,1)*f(:,1) + Gr(:,:,2)*f(:,2);
ucv(3,bxs) = real(u(2));

end


u = G(:,:,1)*f(:,1) + G(:,:,2)*f(:,2);
disp(real(u'))
disp(real(Gh(:,:,1)*f(:,1) + Gh(:,:,2)*f(:,2))')
disp(real(Gr(:,:,1)*f(:,1) + Gr(:,:,2)*f(:,2))')

%% Try ust calculating velocity first?? Try
% u = [0;0;0];
% u2 = u;
% f = [1;0;0];
% f(:,2) = [-1+1;0;0];
% x0 = [.5;.5;.5];
% x02 = [.75;.5;.5];
% x0(:,2) = x02;
% x0(:,3) = [.625;.5;.5];
% V = 1;
% it = 1;
% uSR = zeros(41*41*41-1,1);
% uSM = uSR;
% ur = uSM;
% r = uSR;
% Gij2 = zeros(3,3,2);
% Gij3 = Gij2;
% 
% for i = -25:25
%     for j = -25:25
%         for k = -25:25
%             for pp = 1:1
%             if (i==0 && j==0 && k==0 && pp==1);continue;end
%             if(pp==1);it = it+1;end
% %           Vectors (real and reciprocal)            
%             p = [i;j;k];
%             kv= 2*pi*[i;j;k];
%             
% %           Real part
%             A = shortG(p + x0(:,1) - x0(:,pp));
% %           Reciprocal part
%             B = specG(kv);
%             if (i==0 && j==0 && k==0);B(:)=0;end
%             
% %           Add in contributions
%             tmp =  A*f(:,pp);
%             uSR(it) = tmp(1);
%             u = u + A*f(:,pp);
%              
%             tmp = exp(-norm(kv)^2/4/eps1^2)/V*B*f(:,pp)*exp(-1i*dot(kv,x0(:,1) - x0(:,pp)));
%             uSM(it) = tmp(1);
%             u = u + B/V*f(:,pp)*exp(-1i*dot(kv,x0(:,1) - x0(:,pp)));
%             if isnan(u(1))
%                 disp('?')
%             end
%             
%             r(it) = norm(p);
% %           Normal sum??
%             tmp = Gij(p + x0(:,1) - x0(:,pp))*f;
%             ur(it) = tmp(1);
% %             u2 = u2 + Gij(p + x0(:,1) - x0(:,pp))*f(:,pp);
%             
% %           Just doing the Green's fn
%             Gij2(:,:,pp) = Gij2(:,:,pp) ...
%                          + B/V*exp(-1i*dot(kv,x0(:,1) - x0(:,pp))) ...
%                          + A;
%             Gij3(:,:,pp) = Gij3(:,:,pp) + Gij(p + x0(:,1) - x0(:,pp));
%             end
%         end
%     end
% end
% % Self-interaction
% us = 8*eps1/sqrt(pi)*f(:,1);
% u = u -us;
% uwu = Gij2(:,:,1)*f(:,1) + Gij2(:,:,2)*f(:,2) - us;
% 
% 
% %% Velocity purely with Fourier, e.g. in Poz
% V = 1;
% it = 0;
% uSR = zeros(41*41*41-1,1);
% uSM = uSR;
% ur = uSM;
% r = uSR;
% G = zeros(3,3);
% x0 = [0.5;0.5;0.5];
% it = 0;
% for i = -20:20
%     for j = -20:20
%         for k = -20:20
%             if (i==0 && j==0 && k==0);continue;end
%             it = it+1;
%             kl= 2*pi*[i;j;k];
%             Gijl = fourG(kl)/dot(cross([1;0;0],[0;1;0]),[0;0;1]);
%             G = G + Gijl*exp(-1i*dot(kl,x0));
%             tst(it) = Gijl(1)*exp(-1i*dot(kl,x0));
%             r(it) = norm(kl);
%         end
%     end
% end
% 
% 
% %% Do for the Green's function 2/r
% G = 0;
% G2= 0;
% x0 = [.5;.5;.5];
% V = 1;
% it = 0;
% GSR = zeros(41*41*41-1,1);
% GSM = GSR;
% Gr = GSM;
% r = GSR;
% eps1 = 0.1;
% for i = -20:20
%     for j = -20:20
%         for k = -20:20
%             if (i==0 && j==0 && k==0);continue;end
%             it = it+1;
%             p = [i;j;k];
%             kv= [i;j;k]*2*pi;
%             
%             G = G + rspecG(kv)/V;
%             G = G + rshortG(p);
%             G2 = G2 + 2/norm(p);
%             r(it) = norm(p);
%             GSM(it) = rspecG(kv)/V;
%             GSR(it) = rshortG(p);
%             Gr(it) = 2/norm(p);
%         end
%     end
% end

%% Green's function functions

% Pozrikidis real portion
function G = shortG(x)
global eps1
r = norm(x);
er = r*eps1;
C = erfc(er) + 2/sqrt(pi)*(2*er^2 - 3)*er*exp(-er^2);
D = erfc(er) + 2/sqrt(pi)*(1 - 2*er^2)*er*exp(-er^2);
G = eye(3)*C/r + x*x'*D/r^3;
end

% Poz Fourier portion
function G = specG(k)
global eps1
kn = norm(k);
w = kn/eps1;
G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2 + 1/8)*(kn^2*eye(3) - k*k')*exp(-w^2/4);
end

% Hasimotos
function G = HshortG(x)
global eps1
r = norm(x);
er = r*eps1;
C = erfc(er) - 2/sqrt(pi)*er*exp(-er^2);
D = erfc(er) + 2/sqrt(pi)*er*exp(-er^2);
G = eye(3)*C/r + x*x'*D/r^3;
end

function G = HspecG(k)
global eps1
k2 = k(1)*k(1) + k(2)*k(2) + k(3)*k(3);
c = k2/4/eps1/eps1;
d = 8*pi*(1+c)*exp(-c)/k2/k2;
G = d*(k2*eye(3) - k*k');

% kn = norm(k);
% w = kn/eps1;
% G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2      )*(kn^2*eye(3) - k*k')*exp(-w^2/4);
end



%% Random strange other Green's functions


% function G = rshortG(x)
% global eps1
% r = norm(x);
% er = eps1*r;
% G = 2/r*(2*er/sqrt(pi)*exp(-er^2)*(er^2-2) + erfc(er));
% end

% function G = rspecG(k)
% global eps1
% kn = norm(k);
% w = kn/eps1;
% G = 8*pi/kn^2*(1/w^4 + 1/4/w^2 + 1/8)*exp(-w^2/4)*1; % 1 for now since no sum
% end

% More general version: screens, takes derivatives, and only plugs in the 
% form of the split at the end, useful for spectral - real
function G = shortG_der(x)
global eps1
r = norm(x);
er = r*eps1;
fp = erfc(er) - 2*er*exp(-er^2)/sqrt(pi);
fppr= -4*er*exp(-er^2)/sqrt(pi) + 4*er^3*exp(-er^2)/sqrt(pi);
C = fppr + fp;
D = fp - fppr;
G = eye(3)*C/r + x*x'*D/r^3;
end

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


% function Gh = specG(k,x0) %%% STRANGE STUFF HERE
% global eps1
% kn = norm(k);
% w = kn/eps1;
% G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2 + 1/8)*(kn^2*eye(3) - k*k');
% end

% function Gh = specG(k) %%% STRANGE STUFF HERE
% global eps1
% w = norm(k)/eps1;
% kn = norm(k);
% Gh = 8*pi/(kn^2)*(eye(3) - k*k'/kn^2)*(1 + 1/4*w^2 + 1/8*w^4)*exp(-1/4*w^2);
% end



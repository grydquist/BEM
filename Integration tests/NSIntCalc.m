%% Reference Shape Definition

% Order of the force and velocity transforms
p = 16;

% Constants for phi
np = 2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';

% Constants for theta
nt = (p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);

h = sqrt(4*pi)/nt;

% Weights for rotated integral calculation
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:p
        tmpL = legendre(j,xs(i));
        ws(i) = ws(i) + tmpL(1)/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

q = 2*p;
ntf = (q+1);
npf = 2*(q+1);

dphif = 2*pi/npf;
phif = 0:dphif:dphif*(npf-1)';
[xsf,wgf] = lgwt(ntf,-1,1);
thtf = acos(xsf);
[phf,thf] = meshgrid(phif,thtf);


% Harmonics evaluated at the given thetas/phis
Yt = SpHarmTNew(p,th,ph);

% Spherical harmonic coefficients defining the shape of ref & def:
xmnR = zeros((p+1)^2,1);
xmnR(1) = 2*sqrt(pi);

% Get harmonics for the individual parts of the position vector (would be
% good to do this for the reference too, but don't need to at this point)
rrr = SpHReconst(xmnR,Yt);
x1 = rrr.*sin(th).*cos(ph);
x2 = rrr.*sin(th).*sin(ph);
x3 = rrr.*cos(th);
x1mn = SpT(Yt,x1,tht,phi);
x2mn = SpT(Yt,x2,tht,phi);
x3mn = SpT(Yt,x3,tht,phi);

x1mn(abs(x1mn)<1e-12) = 0;
x2mn(abs(x2mn)<1e-12) = 0;
x3mn(abs(x3mn)<1e-12) = 0;

% Get the coordinates on grid
x = zeros(3,nt,np);
x(1,:,:) = SpHReconst(x1mn,Yt);
x(2,:,:) = SpHReconst(x2mn,Yt);
x(3,:,:) = SpHReconst(x3mn,Yt);

%% Calc the integral
b = [0,0,0];
xcr = [0,0,1.1];
for i = 1:nt
    for j = 1:np
        xi = [x(1,i,j),x(2,i,j),x(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        b = b + [1;1;1]'*Tij(r',n)*dphi*wg(i);
%         b = b + [1;1;1]'*Gij(r')*dphi*wg(i);
%         b = b + [1;1;1]'*Tij(r',n)*dphi*ws(i)*sin(tht(i));
%         b = b + [1;1;1]'*Gij(r')*dphi*ws(i)*sin(tht(i));
    end
end

% bs(p) = b(3);

%% Integral using QBX
% This seems quite close to working

xc = [0,0,2];
xx = [0,0,1];
% Number of coeffs
c = 10;
zs= zeros(3, c + 1);

rx = xx-xc;
% Loop to get the QBX coeffs of order n
for k = 0:c
for i = 1:nt
    for j = 1:np
        xi = [x(1,i,j),x(2,i,j),x(3,i,j)];
        n  = xi;
        r  = xi - xx;
        ry = xi - xc;
        rym = norm(ry);
        rxm = norm(rx);
        
        costh = dot(rx,ry)/(rxm*rym);% to optimize, switch up loop order
        tmpL = legendre(k,costh);
        lP = tmpL(1);
        lD = legendreD(k,costh);
        
        zs(:,k+1) = zs(:,k+1) + ([1,1,1]*(eye(3)*lP/(rym^(k+1)) ...
                  - r'*((k+1)*(-ry)/(rym^(k+3))*lP ... 
                  + lD*(rx/(rxm*rym) + dot(rx,ry)*-ry/(rxm*rym^3)) ...
                  / (rym^(k+1)))))'*dphi*wg(i);
%         zs(:,k+1) = zs(:,k+1) + ([1;1;1]'*Tij(r',n)*norm(r)*dphi*wg(i)*lP/(norm(ry)^(k+1)))';
%         b2 = b2 + [1;1;1]'*Tij(r',n)*dphi*ws(i)*sin(tht(i));
    end
end
end

b2 = [0,0,0];
for k = 0:c
   b2 = b2 + zs(:,k+1)'*rxm^(k);
end

if exist('xax','var')
    xax = 1:max([length(xax),c]);
else
    xax=1;
    xax = 1:max([length(xax),c]);
end
    
bsI(c) = b2(1);
b2(1)

%% Integral with just the interpolation
% Seems unsuitable b/c principal value??
c = 6;
xcs = linspace(1,1.5,c);
bcs = zeros(c,1);

for l = 1:c
    
xcr = [0,0,xcs(l)];
bc = [0,0,0];
for i = 1:nt
    for j = 1:np
        xi = [x(1,i,j),x(2,i,j),x(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        if l>1
%             b = b + [1;1;1]'*Gij(r')*dphi*wg(i);
            bc = bc + [1;1;1]'*Tij(r',n)*dphi*wg(i);
        else
%             b = b + [1;1;1]'*Gij(r')*dphi*ws(i)*sin(tht(i));
            bc = bc + [1;1;1]'*Tij(r',n)*dphi*ws(i)*sin(tht(i));
        end
    end
end
    bcs(l) = bc(1);
    
end

a = polyfit(xcs,bcs,c-1);
bpol = polyval(a,1.1);

%% Method found in Zhao: in mine, I'll rotate for precomp I think

%%%%%%%% THIS IS PROMISING. INSTEAD OF MASK/ PART OF UNITY, JUST DO
% REFINED INTEGRAL AT PATCH (MAYBE REFINE MORE, make bigger)

%%% DOESN'T MAKE SENSE TO OVERLAP INTS AS ARE??? but cutoff does need
%%% thought
% Instead of just saying = 0, just do two separate integrals
t = flip((xs + 1)/2);
th1 = tht(3)*pi/acos(2*t(1)-1);%pi/(p)^(1/2);
thtz = acos(2*(flip(t)-0.5))*th1/pi;
% thtz = flip(th1/2*xs + th1/2); % for no cos trans
msk = exp((2*exp(-1./t)./(t-1)));
msk(:)=1;
% Calculate the masking function for use in the 2nd, nonsingular integ
mskNS = zeros(nt,1);
for i = 1:nt
    if tht(i) < tht(3)-.001%th1 %% this exact matching maybe doesn't male sense
        mskNS(i) = 1;%exp((2*exp(-1./(tht(i)/th1))./((tht(i)/th1)-1)));
    else
        mskNS(i) = 0;
    end
end

% Calculate first integral

I1 = [0,0,0];
xcr = [0,0,1.1];
sm1 = squeeze(x(1,:,:));
for i = 1:nt
    for j = 1:np
        xi = [x(1,i,j),x(2,i,j),x(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        I1 = I1 + [1;1;1]'*Tij(r',n)*dphi*wg(i)*(1-mskNS(i));
%         I1 = I1 + [1;1;1]'*Gij(r')*dphi*wg(i)*(1-mskNS(i));
        tmp =     Tij(r',n)*sin(tht(i))*(1-mskNS(i));%*dphi*wg(i);
        sm1(i,j) = tmp(1,1);
%         b = b + [1;1;1]'*Tij(r',n)*dphi*ws(i)*sin(tht(i));
%         b = b + [1;1;1]'*Gij(r')*dphi*ws(i)*sin(tht(i));
    end
end

% Second integral
x2 = x;
[phz,thz] = meshgrid(phi,thtz);
Yz = SpHarmTNew(p,thz,phz);
x2(1,:,:) = SpHReconst(x1mn,Yz);
x2(2,:,:) = SpHReconst(x2mn,Yz);
x2(3,:,:) = SpHReconst(x3mn,Yz);

I2 = [0,0,0];
for i = 1:nt
    for j = 1:np
        xi = [x2(1,i,j),x2(2,i,j),x2(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        I2 = I2 + [1;1;1]'*Tij(r',n)*dphi*wg(i)*th1/pi*sin(thtz(i))/sin(thtz(i)*pi/th1)*msk(i); 
%         I2 = I2 + [1;1;1]'*Tij(r',n)*dphi*wg(i)*sin(thtz(i))*th1/2*msk(i); % For no cos trans
        
        tmp =              Tij(r',n)*sin(thtz(i))*msk(i);
        sm1(i,j) = tmp(1,1);
    end
end

bz = I1 + I2;

%% Redoing the above so that the integrals make a bit more sense
% Local refinement in the highly varying part for higher acc.
% Uses q = 2*p in patch

% Edge of refined polar patch
rho = pi/(5*p)^(1/2);

% Make the new meshes
thtzf = rho/2*(xsf + 1);%rho/pi*acos(xs);
thtzc = xs*(pi-rho)/2 + (pi + rho)/2;%(pi-rho)/pi*acos(xs) + rho;
[tr,thzf] = meshgrid(phif,thtzf);
[ph,thzc] = meshgrid(phi,thtzc);

Yzf = SpHarmTNew(q,thzf,phf);
Yzc = SpHarmTNew(p,thzc,ph);

xzf = zeros(3,ntf,npf);
xzf(1,:,:) = SpHReconst(x1mn,Yzf,p);
xzf(2,:,:) = SpHReconst(x2mn,Yzf,p);
xzf(3,:,:) = SpHReconst(x3mn,Yzf,p);

xzc = x;
xzc(1,:,:) = SpHReconst(x1mn,Yzc);
xzc(2,:,:) = SpHReconst(x2mn,Yzc);
xzc(3,:,:) = SpHReconst(x3mn,Yzc);


% Something appears to be seriously wrong with the accuracy of these
% integrals???
% When rho = pi/2, should be more accurate, but it's much less accurate
I1 = 0;
I2 = 0;
xcr = [0,0,1.1732];

sm1 = squeeze(x(1,:,:));
sm2 = squeeze(x(1,:,:));
for i = 1:ntf
    for j = 1:npf
        xi = [xzf(1,i,j),xzf(2,i,j),xzf(3,i,j)];
        n  = xi;
        r  = xi - xcr;
%       It's for some reason much more accurate to not do cos transfo?
        I1 = I1 + [1;1;1]'*Tij(r',n)*sin(thtzf(i))*(rho/2)*wgf(i)*dphif;
        
        tmp = Tij(r',n)*sin(thtzf(i));
        sm1(i,j) = tmp(1,1);
    end
end
        
for i = 1:nt
    for j = 1:np
        xi = [xzc(1,i,j),xzc(2,i,j),xzc(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        I2 = I2 + [1;1;1]'*Tij(r',n)*wg(i)*sin(thtzc(i))*(pi-rho)/2*dphi;
        
        tmp = Tij(r',n)*sin(thtzc(i));
        sm2(i,j) = tmp(1,1);
    end
end

disp(I1 + I2)

%% Gij test -- This appears to work at a single point
% xc = [0,0,2];
% xx = [0,0,1.5];
% xi = [x(1,1,1),x(2,1,1),x(3,1,1)];
% r  = xi - xx;
% ry = xi - xc;
% rym = norm(ry);
% rx = xx - xc;
% rxm = norm(rx);
% 
% costh = dot(rx,ry)/(rxm*rym);
% 
% b3 = Gij(r);
% c = 12;
% 
% zs = zeros(3, 3, c + 1);
% zs2= zeros(3, c + 1);
% % Calculating via expansion around xc
% for k = 0:c
%     tmpL = legendre(k,costh);
%     lP = tmpL(1);
%     lD = legendreD(k,costh);
%     zs(:,:,k+1) = (eye(3)*lP/(rym^(k+1)) ...
%                   - r'*((k+1)*(-ry)/(rym^(k+3))*lP ... 
%                   + lD*(rx/(rxm*rym) + dot(rx,ry)*-ry/(rxm*rym^3)) ...
%                   / (rym^(k+1))));
%     zs2(:,k+1) = ((k+1)*(-ry)/(rym^(k+3))*lP ... 
%                   + lD*(rx/(rxm*rym) + dot(rx,ry)*-ry/(rxm*rym^3)) ...
%                   / (rym^(k+1)));
% end 
% 
% b4 = b3;
% b4(:,:) = 0;
% b5 = [0,0,0];
% 
% for k = 0:c
%    b4 = b4 + zs(:,:,k+1)'*rxm^k;
%    b5 = b5 + zs2(:,k+1)'*rxm^k;
% end
% 
% 
% 






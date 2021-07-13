%% Reference Shape Definition
%clustered seems to be better, as it is better at resolving rapidly varying
%part (also in Zhao they mention cluster)
% Order of the force and velocity transforms
p = 12;
q = 2*p;

% Constants for phi
np = 2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';

% Constants for theta
nt = (p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);

myfacs = zeros(2*q+1,1);
for i = 0:2*q
    myfacs(i+1) = factorial(i);
end    

% Weights for rotated integral calculation
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:p
        tmpL = legendre(j,xs(i));
        ws(i) = ws(i) + tmpL(1)/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end


ntf = (q+1);
npf = 2*(q+1);

dphif = 2*pi/npf;
phif = 0:dphif:dphif*(npf-1)';
[xsf,wgf] = lgwt(ntf,-1,1);
thtf = acos(xsf);
[phf,thf] = meshgrid(phif,thtf);

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmTNew(p,th,ph);
Ytf = SpHarmTNew(q,thf,phf);

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

% [x1mn,x2mn,x3mn] = RBCcoeffs(Yt,th,ph);

% Get the coordinates on grid
x = zeros(3,nt,np);
x(1,:,:) = SpHReconst(x1mn,Yt);
x(2,:,:) = SpHReconst(x2mn,Yt);
x(3,:,:) = SpHReconst(x3mn,Yt);
xf = zeros(3,ntf,npf);
xf(1,:,:) = SpHReconst(x1mn,Ytf,p);
xf(2,:,:) = SpHReconst(x2mn,Ytf,p);
xf(3,:,:) = SpHReconst(x3mn,Ytf,p);

%% Find the minimum distance, rotate around this point
% xcr = [0,0,1.01];
xcr = [1,1,1]/norm([1,1,1]) + .01;
% xcr = [.81,.81,.42]+.05;

% shortcut, will need to add in later
% xrot = xcr/norm(xcr);
% phir = atan2(xrot(2),xrot(1));
% thtr = acos(xrot(3));


% Instead of finding exact point, just pick the closest point...
% Maybe try using xf instead as well
h = 1;
for i = 1:ntf
    for j = 1:npf
        d = norm(xf(:,i,j) - xcr');
        h = min([h,d]);
        if h == d 
            indi = i;
            indj = j;
        end
    end
end
phir = phif(indj);
thtr = thtf(indi);
disp(h)

x1mnr = SpHRot(x1mn,phir,-thtr,-phir, myfacs);
x2mnr = SpHRot(x2mn,phir,-thtr,-phir, myfacs);
x3mnr = SpHRot(x3mn,phir,-thtr,-phir, myfacs);

xcf = zeros(3,nt,np);
xcf(1,:,:) = real(SpHReconst(x1mnr,Yt));
xcf(2,:,:) = real(SpHReconst(x2mnr,Yt));
xcf(3,:,:) = real(SpHReconst(x3mnr,Yt));

%% Calculate the integral centered on this point

rho =pi/6;% pi/(4*p)^(1/2);
h = sqrt(4*pi)/nt;

% Make the new meshes
thtzf = rho/2*(xsf + 1);%rho/pi*acos(xs);
thtzc = xs*(pi-rho)/2 + (pi + rho)/2;%(pi-rho)/pi*acos(xs) + rho;
[tr,thzf] = meshgrid(phif,thtzf);
[ph,thzc] = meshgrid(phi,thtzc);

Yzf = SpHarmTNew(q,thzf,phf);
Yzc = SpHarmTNew(p,thzc,ph);

xzf = zeros(3,ntf,npf);
xzf(1,:,:) = real(SpHReconst(x1mnr,Yzf,p));
xzf(2,:,:) = real(SpHReconst(x2mnr,Yzf,p));
xzf(3,:,:) = real(SpHReconst(x3mnr,Yzf,p));

xzc = zeros(3,nt,np);
xzc(1,:,:) = real(SpHReconst(x1mnr,Yzc,p));
xzc(2,:,:) = real(SpHReconst(x2mnr,Yzc,p));
xzc(3,:,:) = real(SpHReconst(x3mnr,Yzc,p));


I1 = 0;
I2 = 0;

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

bz = I1 + I2;
a(p) = bz(1);

%% Normal integral
b = [0,0,0];
for i = 1:ntf
    for j = 1:npf
        xi = [xf(1,i,j),xf(2,i,j),xf(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        b = b + [1;1;1]'*Tij(r',n)*wgf(i)*sin(thtf(i))*dphif;
        
    end
end

%% Sinh transform

% Not grid spacing, but separation distance
h = sqrt(pi)/10/nt;
k = -0.5*log(rho/h + sqrt((rho/h)^2+1));
z = xsf;
thth = h*sinh(k*(z-1));

[phf,thh] = meshgrid(phif,thth);

Yh = SpHarmTNew(q,thh,phf);

xh = zeros(3,ntf,npf);
xh(1,:,:) = real(SpHReconst(x1mnr,Yh,p));
xh(2,:,:) = real(SpHReconst(x2mnr,Yh,p));
xh(3,:,:) = real(SpHReconst(x3mnr,Yh,p));

I1 = 0;
for i = 1:ntf
    for j = 1:npf
        xi = [xh(1,i,j),xh(2,i,j),xh(3,i,j)];
        n  = xi;
        r  = xi - xcr;
        I1 = I1 + [1;1;1]'*Tij(r',n)*sin(thth(i))*h*-k*cosh(k*xsf(i)-k)*wgf(i)*dphif;
        
        tmp = Tij(r',n)*sin(thth(i));
        sm1(i,j) = tmp(1,1);
    end
end
        
I2 = 0;
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

disp(I1+I2)

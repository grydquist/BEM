% Differences: get the deformed state info first. With that, calculate the
% velocity. Only need to get material constants once.

%% Material Constants

%notes: lam=1 corresponds to A = 4pi identity which falls out from the
%governing eq.

% Viscosity outside
mu = 1;
% Viscosity ratio
lam = 1;%/5;

% Deformation resistance constants
% B = 12.4;
% C = 200;
B = .005;
C = 100;

% Bending modulus
Eb = 0.0669;

% Total time steps
NT = 1400;
dt = 0.005;

% Velocity and gradient
U = [0;0;0];

shft = [cos(pi/4),0,sin(pi/4);0,1,0;-sin(pi/4),0,cos(pi/4)];
dU = [0,0,1;0,0,0;.0,0,0];
% dU = shft'*[0,0,1;0,0,0;1,0,0]*shft;%/4/pi*shft;
%% Reference Shape Definition

% Order of the force and velocity transforms
p = 8;

% To de-alias the force, we need to get a finer grid. The factor is fali
fali = 2;

% Order of SpH for the fine grid
q = p*fali;

% I call factorial a lot and it makes the most sense to do all the ones
% I'll need right off the bat and store them (goes up to 2p)
myfacs = zeros(2*q+1,1);
for i = 0:2*q
    myfacs(i+1) = factorial(i);
end    

% How many coefficients per dimension?
ftot = (p+1)^2;

% Number of points times more than needed for the integrals
gf = 1;

% Constants for phi
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';

% Constants for theta
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);

% Making the fine grid
ntf = gf*(q+1);
npf = gf*2*(q+1);

dphif = 2*pi/npf;
phif = 0:dphif:dphif*(npf-1)';
[xsf,wgf] = lgwt(ntf,-1,1);
thtf = acos(xsf);
[phf,thf] = meshgrid(phif,thtf);

% Weights for rotated integral calculation
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:gf*p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmT(p,th,ph,myfacs);

% Harmonics evaluated at fine grid
Ytf = SpHarmT(q,thf,phf,myfacs);

% Spherical harmonic coefficients defining the shape of ref & def:
xmnR = zeros((p+1)^2,1);
xmnR(1) = 2*sqrt(pi);

% Get harmonics for the individual parts of the position vector (would be
% good to do this for the reference too, but don't need to at this point)
rrr = SpHReconst(xmnR,Yt);
x1 = rrr.*sin(th).*cos(ph);
x2 = rrr.*sin(th).*sin(ph);
x3 = rrr.*cos(th);

% For the RBC, fully defined at p = 5
[x1mn,x2mn,x3mn] = RBCcoeffs(Yt,th,ph);
xmns = zeros((p+1)^2,3,NT);
xmns(:,1,1) = x1mn;
xmns(:,2,1) = x2mn;
xmns(:,3,1) = x3mn;

x1mn(abs(x1mn)<1e-12) = 0;
x2mn(abs(x2mn)<1e-12) = 0;
x3mn(abs(x3mn)<1e-12) = 0;

% x1mn = SpT(Yt,x1,th,ph);
% x2mn = SpT(Yt,x2,th,ph);
% x3mn = SpT(Yt,x3,th,ph);

ix3mn = x3mn;
J = zeros(ntf,npf);

% Get the coordinates on grid
x = zeros(3,nt,np);
x(1,:,:) = SpHReconst(x1mn,Yt);
x(2,:,:) = SpHReconst(x2mn,Yt);
x(3,:,:) = SpHReconst(x3mn,Yt);

%% Preallocation of all the solid things I need

% Fine grid
xf = zeros(3,ntf,npf);

% Total derivatives, in Cartesian
dxp = zeros(3,ntf,npf);
dxt = dxp;

% Second order
dxp2= dxp;
dxt2= dxp;
dxtp= dxp;

% Third order
dxp3 = dxp;
dxt3 = dxp;
dxt2p= dxp;
dxtp2= dxp;

% Fourth order
dxp4 = dxp;
dxtp3 = dxp;
dxt2p2 = dxp;
dxt3p = dxp;
dxt4 = dxp;

% Normal Vector (inward)
nk = zeros(3,ntf,npf);

% Moment
dmabt = zeros(2,2);
dmabp = zeros(2,2);
dmab2 = zeros(2,2);

dtau  = zeros(2,2);
dq = zeros(2,1);
bm = dtau;
bn = bm;

myf = zeros(3,ntf,npf);
Nmyf = zeros(3,nt,np);
fab = myf;

% Energy spectrums

Est = zeros(1,q+1);
Es1 = zeros(1,q+1);
Es2 = zeros(1,q+1);
Es3 = zeros(1,q+1);
Ex = zeros(1,p+1);
Eu = Ex;

c1R = zeros(3,ntf,npf);
c2R = c1R;
c1tR= c1R;
c1pR= c1R;
c2tR= c1R;
c2pR= c1R;
kR = zeros(ntf,npf);
kdR = zeros(2,ntf,npf);
kd2R = zeros(3,ntf,npf);

% Get all the derivatives of the Spherical Harmonics up front
Ytd1 = zeros(ntf,npf,(p+1)^2);
Ytd2 = Ytd1;
Ytd3 = Ytd1;
Ytd4 = Ytd1;

% Also needed on the coarser grid
Ytd1_c = zeros(nt,np,(p+1)^2);

it = 0;
im = 0;
for n = 0:p
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
        
%       Theta derivative
        Ytd1(:,:,it) = ThetDer(Ytf,phf,n,m,1);
        Ytd1_c(:,:,it)= ThetDer(Yt,ph,n,m,1);
        
%       Second theta derivative 
        Ytd2(:,:,it) = ThetDer(Ytf,phf,n,m,2);

%       Third theta derivative
        Ytd3(:,:,it) = ThetDer(Ytf,phf,n,m,3);
        
%       Fourth theta derivative
        Ytd4(:,:,it) = ThetDer(Ytf,phf,n,m,4);
    end
end

%% Fluid Preallocation

A = zeros(3*nt*np, 3*ftot);
At = zeros(3,3);
v = At;
b = zeros(3*nt*np,1);
bt = zeros(3,1);
thet = zeros(nt,np);
phit = thet;
A2 = zeros(3*ftot,3*ftot);
b2 = zeros(3*ftot,1);
nkg = zeros(3,nt,np);
Jg  = zeros(nt,np);
dxtg = zeros(3,nt,np);
dxpg = dxtg;
xcg = zeros(3,nt,np);
ua = zeros(3,nt,np);
ua11 = ua;
    
%% GIF writing
h = figure(1);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'pull.gif';

% Output thing I don't want to run every timestep
tmpt = linspace(0,pi,100);tmpp = linspace(0,pi,100);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmT(p,tmpth,tmpph,myfacs);

trcpnt = zeros(1,NT);
trcpntz = zeros(1,NT);
trcvel = zeros(1,NT);
trcvelz = zeros(1,NT);

% Spherical harmonic evaluated at right hand side of sphere
Ytrc = SpHarmT(p,pi/2,0,myfacs);

%% Time stepping loop
for cts = 1:NT
%% Calculation of residual force at convected points
% After calculation, do SpHarm transform
% if(cts==25);dU(:)=0;end %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Get Geometric information for current configuration (Perhaps the best way
% to do this is with a subroutine that takes in the SpHarm consts for the
% def, and infor for the ref, and material constants and just spits out the
% forces. In FORTRAN this information would be packaged up in objects, and
% we could just pass pointers to the objects!)

ih = 0;
it = 0;
xf(:) = 0;
dxt(:) = 0;
dxp(:) = 0;
dxt2(:) = 0;
dxp2(:) = 0;
dxtp(:) = 0;
dxt3(:) = 0;
dxp3(:) = 0;
dxt2p(:) = 0;
dxtp2(:) = 0;
dxp4(:) = 0;
dxtp3(:) = 0;
dxt2p2(:) = 0;
dxt3p(:) = 0;
dxt4(:) = 0;

% Calcaulating surface derivatives
% Loop over harmonics
for n = 0:p
    ih = ih+1;
%   This is assuming that the deformed and reference are at the same tht
%   and phi, which they should be, as they aren't angles necessarily
    Ypcur = Ytf{ih};
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
%       Current value of harmonic coefficient        
        f1 = x1mn(it);
        f2 = x2mn(it);
        f3 = x3mn(it);
%       Current values of harmonic at degree and order
        Ym = squeeze(Ypcur(im,:,:));
%       Theta derivative
        td1 = Ytd1(:,:,it);
        
%       Distance from zero
        xf(1,:,:) = squeeze(xf(1,:,:)) + f1*Ym;
        xf(2,:,:) = squeeze(xf(2,:,:)) + f2*Ym;
        xf(3,:,:) = squeeze(xf(3,:,:)) + f3*Ym;
        xf = real(xf);
        
%       Get the derivative in theta and add in
        dxt(1,:,:) = squeeze(dxt(1,:,:)) + f1*td1;
        dxt(2,:,:) = squeeze(dxt(2,:,:)) + f2*td1;
        dxt(3,:,:) = squeeze(dxt(3,:,:)) + f3*td1;
        
%       Get the derivative in phi and add in
        dxp(1,:,:) = squeeze(dxp(1,:,:)) + f1*1i*m*Ym;
        dxp(2,:,:) = squeeze(dxp(2,:,:)) + f2*1i*m*Ym;
        dxp(3,:,:) = squeeze(dxp(3,:,:)) + f3*1i*m*Ym;
        
%       Second order derivatives (for curvature
%       Second theta derivative 
        td2 = Ytd2(:,:,it);

%       Second phi derivative
        dxp2(1,:,:) = squeeze(dxp2(1,:,:)) + f1*-m^2*Ym;
        dxp2(2,:,:) = squeeze(dxp2(2,:,:)) + f2*-m^2*Ym;
        dxp2(3,:,:) = squeeze(dxp2(3,:,:)) + f3*-m^2*Ym;

%       Mixed derivative
        dxtp(1,:,:) = squeeze(dxtp(1,:,:)) + f1*1i*m*td1;
        dxtp(2,:,:) = squeeze(dxtp(2,:,:)) + f2*1i*m*td1;
        dxtp(3,:,:) = squeeze(dxtp(3,:,:)) + f3*1i*m*td1;

%       Second theta derivativeS
        dxt2(1,:,:) = squeeze(dxt2(1,:,:)) + f1*td2;
        dxt2(2,:,:) = squeeze(dxt2(2,:,:)) + f2*td2;
        dxt2(3,:,:) = squeeze(dxt2(3,:,:)) + f3*td2;

%       Third order derivatives (for curvature derivatives)
%       Third phi derivative
        dxp3(1,:,:) = squeeze(dxp3(1,:,:)) + f1*-1i*m^3*Ym;
        dxp3(2,:,:) = squeeze(dxp3(2,:,:)) + f2*-1i*m^3*Ym;
        dxp3(3,:,:) = squeeze(dxp3(3,:,:)) + f3*-1i*m^3*Ym;

%       Mixed derivatives
        dxt2p(1,:,:) = squeeze(dxt2p(1,:,:)) + f1*1i*m*td2;
        dxt2p(2,:,:) = squeeze(dxt2p(2,:,:)) + f2*1i*m*td2;
        dxt2p(3,:,:) = squeeze(dxt2p(3,:,:)) + f3*1i*m*td2;
        
        dxtp2(1,:,:) = squeeze(dxtp2(1,:,:)) + f1*-m^2*td1;
        dxtp2(2,:,:) = squeeze(dxtp2(2,:,:)) + f2*-m^2*td1;
        dxtp2(3,:,:) = squeeze(dxtp2(3,:,:)) + f3*-m^2*td1;

%       Third theta derivative
        td3 = Ytd3(:,:,it);
        dxt3(1,:,:) = squeeze(dxt3(1,:,:)) + f1*td3;
        dxt3(2,:,:) = squeeze(dxt3(2,:,:)) + f2*td3;
        dxt3(3,:,:) = squeeze(dxt3(3,:,:)) + f3*td3;
        
        
%       Fourth order derivatives (for shear force derivatives)
%       Fourth phi derivative
        dxp4(1,:,:) = squeeze(dxp4(1,:,:)) + f1*m^4*Ym;
        dxp4(2,:,:) = squeeze(dxp4(2,:,:)) + f2*m^4*Ym;
        dxp4(3,:,:) = squeeze(dxp4(3,:,:)) + f3*m^4*Ym;

%       Mixed derivatives
        dxtp3(1,:,:) = squeeze(dxtp3(1,:,:)) + f1*-1i*m^3*td1;
        dxtp3(2,:,:) = squeeze(dxtp3(2,:,:)) + f2*-1i*m^3*td1;
        dxtp3(3,:,:) = squeeze(dxtp3(3,:,:)) + f3*-1i*m^3*td1;
        
        dxt2p2(1,:,:) = squeeze(dxt2p2(1,:,:)) + f1*-m^2*td2;
        dxt2p2(2,:,:) = squeeze(dxt2p2(2,:,:)) + f2*-m^2*td2;
        dxt2p2(3,:,:) = squeeze(dxt2p2(3,:,:)) + f3*-m^2*td2;
        
        dxt3p(1,:,:) = squeeze(dxt3p(1,:,:)) + f1*1i*m*td3;
        dxt3p(2,:,:) = squeeze(dxt3p(2,:,:)) + f2*1i*m*td3;
        dxt3p(3,:,:) = squeeze(dxt3p(3,:,:)) + f3*1i*m*td3;

%       Fourth theta derivative
        td4 = Ytd4(:,:,it);
        dxt4(1,:,:) = squeeze(dxt4(1,:,:)) + f1*td4;
        dxt4(2,:,:) = squeeze(dxt4(2,:,:)) + f2*td4;
        dxt4(3,:,:) = squeeze(dxt4(3,:,:)) + f3*td4;
    end
end
% These shouldn't be imaginary, but some creeps in        
dxt   = real(dxt);
dxt2  = real(dxt2);
dxt3  = real(dxt3);

dxp   = real(dxp);
dxp2  = real(dxp2);
dxp3  = real(dxp3);

dxtp  = real(dxtp);
dxt2p = real(dxt2p);
dxtp2 = real(dxtp2);

dxp4 = real(dxp4);
dxtp3 = real(dxtp3);
dxt2p2 = real(dxt2p2);
dxt3p = real(dxt3p);
dxt4 = real(dxt4);


for i = 1:ntf
    for j = 1:npf
        
%       Normal vector (inward)
        nk(:,i,j) = cross(dxt(:,i,j),dxp(:,i,j));
        nk(:,i,j) = -nk(:,i,j)./norm(nk(:,i,j));
        nk(:,i,j) = real(nk(:,i,j));
        
%       Jacobian via fundamental forms
        E = dot(dxt(:,i,j),dxt(:,i,j));
        F = dot(dxt(:,i,j),dxp(:,i,j));
        G = dot(dxp(:,i,j),dxp(:,i,j));
        J(i,j) = sqrt(E.*G - F.^2);

%       Calculate mean curvature with some fundamental forms (appendices in
%       Veeranpani paper)
        L = dot(dxt2(:,i,j),-nk(:,i,j));
        M = dot(dxtp(:,i,j),-nk(:,i,j));
        N = dot(dxp2(:,i,j),-nk(:,i,j));
%       Curvature numerator
        D = E*N - 2*F*M + G*L;
        k = 0.5*D/J(i,j)^2;
        
%       Perhaps simpler is that k is the mean of the principle curvatures.

%       First let's do all the cross products we need
        ct2_p = cross(dxt2(:,i,j),dxp(:,i,j));
        ct_tp = cross(dxt (:,i,j),dxtp(:,i,j));
        ctp_p = cross(dxtp(:,i,j),dxp(:,i,j));
        ct_p2 = cross(dxt (:,i,j),dxp2(:,i,j));
        
        ct3_p = cross(dxt3(:,i,j),dxp(:,i,j));
        ct2_tp = cross(dxt2(:,i,j),dxtp(:,i,j));
        ct_t2p = cross(dxt(:,i,j),dxt2p(:,i,j));
        
        ct2p_p = cross(dxt2p(:,i,j),dxp(:,i,j));
        ct2_p2 = cross(dxt2(:,i,j),dxp2(:,i,j));
        ct_tp2 = cross(dxt(:,i,j),dxtp2(:,i,j));
        
        ctp2_p = cross(dxtp2(:,i,j),dxp(:,i,j));
        ctp_p2 = cross(dxtp(:,i,j),dxp2(:,i,j));
        ct_p3 = cross( dxt(:,i,j),dxp3(:,i,j));
        
%       Jacobian Partials        
        Jt = dot(-nk(:,i,j),ct2_p + ct_tp);
        Jp = dot(-nk(:,i,j),ctp_p + ct_p2);

%       Normal vector partials outward
        dnt = (1/J(i,j))*(ct2_p + ct_tp - Jt*-nk(:,i,j));
        dnp = (1/J(i,j))*(ctp_p + ct_p2 - Jp*-nk(:,i,j));
        
%       Second Jacobian derivs
        Jt2 = dot(dnt,ct2_p + ct_tp) ...
            + dot(-nk(:,i,j), ct3_p + 2*ct2_tp + ct_t2p);
        Jtp = dot(dnp,ct2_p + ct_tp) ...
            + dot(-nk(:,i,j), ct2p_p +  ct2_p2 + ct_tp2);
        Jp2 = dot(dnp,ctp_p + ct_p2) ...
            + dot(-nk(:,i,j), ctp2_p + 2*ctp_p2 + ct_p3);
        
%       Second normal derivatives        
        dnt2 = (1/J(i,j))*(-2*Jt*dnt - -nk(:,i,j)*Jt2 + ct3_p + 2*ct2_tp + ct_t2p);
        dnp2 = (1/J(i,j))*(-2*Jp*dnp - -nk(:,i,j)*Jp2 + ct_p3 + 2*ctp_p2 + ctp2_p);
        dntp = (1/J(i,j))*(-Jp*dnt - Jt*dnp - -nk(:,i,j)*Jtp + ct2p_p + ct2_p2 + ct_tp2);
     
%       Fundamental form partials        
        Et  = 2*dot(dxt2(:,i,j),dxt(:,i,j));
        Ep  = 2*dot(dxtp(:,i,j),dxt(:,i,j));
        Et2 = 2*(dot(dxt3(:,i,j),dxt(:,i,j)) + dot(dxt2(:,i,j),dxt2(:,i,j)));
        Ep2 = 2*(dot(dxtp2(:,i,j),dxt(:,i,j))+ dot(dxtp(:,i,j),dxtp(:,i,j)));
        Etp = 2*(dot(dxt2p(:,i,j),dxt(:,i,j))+ dot(dxt2(:,i,j),dxtp(:,i,j)));
        
        Gt  = 2*dot(dxtp(:,i,j),dxp(:,i,j));
        Gp  = 2*dot(dxp (:,i,j),dxp2(:,i,j));
        Gt2 = 2*(dot(dxt2p(:,i,j),dxp(:,i,j))+ dot(dxtp(:,i,j),dxtp(:,i,j)));
        Gp2 = 2*(dot(dxp3(:,i,j),dxp(:,i,j)) + dot(dxp2(:,i,j),dxp2(:,i,j)));
        Gtp = 2*(dot(dxtp2(:,i,j),dxp(:,i,j))+ dot(dxtp(:,i,j),dxp2(:,i,j)));
        
        Ft  = dot(dxt2(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxtp(:,i,j));
        Fp  = dot(dxtp(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxp2(:,i,j));
        Ft2 = dot(dxt3(:,i,j),dxp(:,i,j)) + 2*dot(dxt2(:,i,j),dxtp(:,i,j)) + dot(dxt(:,i,j),dxt2p(:,i,j));
        Fp2 = dot(dxtp2(:,i,j),dxp(:,i,j))+ 2*dot(dxtp(:,i,j),dxp2(:,i,j)) + dot(dxt(:,i,j),dxp3(:,i,j));
        Ftp = dot(dxt2p(:,i,j),dxp(:,i,j))+ dot(dxtp(:,i,j),dxtp(:,i,j)) + dot(dxt2(:,i,j),dxp2(:,i,j)) + dot(dxt(:,i,j),dxtp2(:,i,j));
        
        Lt  = dot(dxt3(:,i,j),-nk(:,i,j))  + dot(dxt2(:,i,j),dnt);
        Lp  = dot(dxt2p(:,i,j),-nk(:,i,j)) + dot(dxt2(:,i,j),dnp);
        Lt2 = dot(dxt4(:,i,j),-nk(:,i,j))  + 2*dot(dxt3(:,i,j),dnt) + dot(dxt2(:,i,j),dnt2);
        Lp2 = dot(dxt2p2(:,i,j),-nk(:,i,j))+ 2*dot(dxt2p(:,i,j),dnp) + dot(dxt2(:,i,j),dnp2);
        Ltp = dot(dxt3p(:,i,j),-nk(:,i,j)) + dot(dxt2p(:,i,j),dnt) + dot(dxt3(:,i,j),dnp) + dot(dxt2(:,i,j),dntp);
        
        Nt  = dot(dxtp2(:,i,j),-nk(:,i,j)) + dot(dxp2(:,i,j),dnt);
        Np  = dot(dxp3(:,i,j),-nk(:,i,j))  + dot(dxp2(:,i,j),dnp);
        Nt2 = dot(dxt2p2(:,i,j),-nk(:,i,j))+ 2*dot(dxtp2(:,i,j),dnt) + dot(dxp2(:,i,j),dnt2);
        Np2 = dot(dxp4(:,i,j),-nk(:,i,j))  + 2*dot(dxp3(:,i,j),dnp) + dot(dxp2(:,i,j),dnp2);
        Ntp = dot(dxtp3(:,i,j),-nk(:,i,j)) + dot(dxtp2(:,i,j),dnp) + dot(dxp3(:,i,j),dnt) + dot(dxp2(:,i,j),dntp);
        
        Mt  = dot(dxt2p(:,i,j),-nk(:,i,j)) + dot(dxtp(:,i,j),dnt);
        Mp  = dot(dxtp2(:,i,j),-nk(:,i,j)) + dot(dxtp(:,i,j),dnp);
        Mt2 = dot(dxt3p(:,i,j),-nk(:,i,j)) + 2*dot(dxt2p(:,i,j),dnt) + dot(dxtp(:,i,j),dnt2);
        Mp2 = dot(dxtp3(:,i,j),-nk(:,i,j)) + 2*dot(dxtp2(:,i,j),dnp) + dot(dxtp(:,i,j),dnp2);
        Mtp = dot(dxt2p2(:,i,j),-nk(:,i,j))+ dot(dxt2p(:,i,j),dnp) + dot(dxtp2(:,i,j),dnt) + dot(dxtp(:,i,j),dntp);

%       Curvature derivatives
%       Numerator derivatives
        Dt = Et*N + E*Nt - 2*Ft*M - 2*F*Mt + Gt*L + G*Lt;
        Dp = Ep*N + E*Np - 2*Fp*M - 2*F*Mp + Gp*L + G*Lp;
        
        Dt2 = Et2*N + 2*Et*Nt + E*Nt2 - 2*Ft2*M - 4*Ft*Mt - 2*F*Mt2 + Gt2*L + 2*Gt*Lt + G*Lt2;
        Dp2 = Ep2*N + 2*Ep*Np + E*Np2 - 2*Fp2*M - 4*Fp*Mp - 2*F*Mp2 + Gp2*L + 2*Gp*Lp + G*Lp2;
        Dtp = Etp*N + Et*Np + Ep*Nt + E*Ntp - 2*Ftp*M - 2*Ft*Mp - 2*Fp*Mt - 2*F*Mtp + Gtp*L + Gt*Lp + Gp*Lt + G*Ltp;
        
        kt = 0.5*Dt*J(i,j)^-2 - 2*k*Jt/J(i,j);
        kp = 0.5*Dp*J(i,j)^-2 - 2*k*Jp/J(i,j);
            
        kt2 = 0.5*Dt2*J(i,j)^-2 - 2*Jt*J(i,j)^-3*Dt +3*J(i,j)^-4*Jt^2*D - J(i,j)^-3*Jt2*D;
        kp2 = 0.5*Dp2*J(i,j)^-2 - 2*Jp*J(i,j)^-3*Dp +3*J(i,j)^-4*Jp^2*D - J(i,j)^-3*Jp2*D;
        ktp = 0.5*Dtp*J(i,j)^-2 - Jt*J(i,j)^-3*Dp - Jp*J(i,j)^-3*Dt +3*J(i,j)^-4*Jt*Jp*D - J(i,j)^-3*Jtp*D;
        
%       Metric tensor
        g = [E,F;F,G];
%       Contravariant g
        gn = inv(g);
%       Partials of COTNRAVARIANT metric tensor (from inverse properties)
        dgt = -gn*[Et,Ft;Ft,Gt]*gn;
        dgp = -gn*[Ep,Fp;Fp,Gp]*gn;
        
        dgt2 = -gn*[Et2,Ft2;Ft2,Gt2]*gn - 2*dgt*[Et,Ft;Ft,Gt]*gn;
        dgp2 = -gn*[Ep2,Fp2;Fp2,Gp2]*gn - 2*dgp*[Ep,Fp;Fp,Gp]*gn;
        dgtp = -gn*[Etp,Ftp;Ftp,Gtp]*gn - 2*dgt*[Ep,Fp;Fp,Gp]*gn;
            
%       Shear transverse
%       Contravariant basis vectors, raise index 
        c1 = dxt(:,i,j)*gn(1,1) + dxp(:,i,j)*gn(1,2);
        c2 = dxt(:,i,j)*gn(2,1) + dxp(:,i,j)*gn(2,2);
        
%       For derivatives of tau/q, we need dc/dtht and dc/dphi
%       Chain rule across raise index operation
        c1t = dxt2(:,i,j)*gn(1,1) + dxtp(:,i,j)*gn(2,1)...
            + dxt(:,i,j)*dgt(1,1) + dxp(:,i,j)*dgt(2,1);
        c1p = dxtp(:,i,j)*gn(1,1) + dxp2(:,i,j)*gn(2,1)...
            + dxt(:,i,j)*dgp(1,1) + dxp(:,i,j)*dgp(2,1);
        c2t = dxt2(:,i,j)*gn(1,2) + dxtp(:,i,j)*gn(2,2)...
            + dxt(:,i,j)*dgt(1,2) + dxp(:,i,j)*dgt(2,2);
        c2p = dxtp(:,i,j)*gn(1,2) + dxp2(:,i,j)*gn(2,2)...
            + dxt(:,i,j)*dgp(1,2) + dxp(:,i,j)*dgp(2,2);
        
%       !!!!!!!!!! HACKY FIRST TIME STEP
        if(cts == 1)
            c1R(:,i,j) = c1;
            c2R(:,i,j) = c2;
            kR(i,j) = k;
            kdR(1,i,j) = kt;
            kdR(2,i,j) = kp;
            kd2R(1,i,j)= kt2;
            kd2R(2,i,j)= kp2;
            kd2R(3,i,j)= ktp;
            
            c1tR(:,i,j) = c1t;
            c1pR(:,i,j) = c1p;
            c2tR(:,i,j) = c2t;
            c2pR(:,i,j) = c2p;
        end
        
%       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
        c111 = dot(dxt2(:,i,j),c1);
        c112 = dot(dxtp(:,i,j),c1);
        c122 = dot(dxp2(:,i,j),c1);
        c222 = dot(dxp2(:,i,j),c2);
        c221 = dot(dxtp(:,i,j),c2);
        c211 = dot(dxt2(:,i,j),c2);
        
%       Christoffel derivatives: need 5 each
        c111t = dot(dxt3(:,i,j) ,c1) + dot(dxt2(:,i,j),c1t);
        c112t = dot(dxt2p(:,i,j),c1) + dot(dxtp(:,i,j),c1t);
        c122t = dot(dxtp2(:,i,j),c1) + dot(dxp2(:,i,j),c1t);
        c221t = dot(dxt2p(:,i,j),c2) + dot(dxtp(:,i,j),c2t);
        c222t = dot(dxtp2(:,i,j),c2) + dot(dxp2(:,i,j),c2t);
        
        c111p = dot(dxt2p(:,i,j),c1) + dot(dxt2(:,i,j),c1p);
        c112p = dot(dxtp2(:,i,j),c1) + dot(dxtp(:,i,j),c1p);
        c211p = dot(dxt2p(:,i,j),c2) + dot(dxt2(:,i,j),c2p);
        c221p = dot(dxtp2(:,i,j),c2) + dot(dxtp(:,i,j),c2p);
        c222p = dot(dxp3(:,i,j) ,c2) + dot(dxp2(:,i,j),c2p);
        
%       Moments and their needed partials
        mab = (k - kR(i,j))*gn;
        dmabt = (kt - kdR(1,i,j))*gn + (k - kR(i,j))*dgt;
        dmabp = (kp - kdR(2,i,j))*gn + (k - kR(i,j))*dgp;
        
        dmab2(1,1) = (kt2 - kd2R(1,i,j))*gn(1,1) + 2*(kt - kdR(1,i,j))*dgt(1,1) + (k - kR(i,j))*dgt2(1,1);
        dmab2(1,2) = (ktp - kd2R(3,i,j))*gn(1,2) + (kt - kdR(1,i,j))*dgp(1,2) ...
                   + (kp  - kdR(2,i,j))*dgt(1,2) + (k  - kR(i,j))*dgtp(1,2);
        dmab2(2,1) = (ktp - kd2R(3,i,j))*gn(2,1) + (kt - kdR(1,i,j))*dgp(2,1) ...
                   + (kp  - kdR(2,i,j))*dgt(2,1) + (k  - kR(i,j))*dgtp(2,1);
        dmab2(2,2) = (kp2 - kd2R(2,i,j))*gn(2,2) + 2*(kp - kdR(2,i,j))*dgp(2,2) + (k - kR(i,j))*dgp2(2,2);
        
%       Transverse Shears (Covariant divergence of moment tensor)
        q1 = Eb*(dmabt(1,1) + dmabp(2,1) + mab(1,1)*(2*c111 + c221) ...
                + mab(2,1)*(2*c112 + c222) + mab(1,2)*c112 + mab(2,2)*c122);
        q2 = Eb*(dmabt(1,2) + dmabp(2,2) + mab(1,2)*(c111 + 2*c221) ...
                + mab(2,2)*(c112 + 2*c222) + mab(1,1)*c211 + mab(2,1)*c221);
            
%       Partials of q, checked seem ok
        dq(1) = Eb*(dmab2(1,1) + dmab2(2,1) ...
              + dmabt(1,1)*(2*c111 + c221) + mab(1,1)*(2*c111t + c221t) ...
              + dmabt(2,1)*(2*c112 + c222) + mab(2,1)*(2*c112t + c222t) ...
              + dmabt(1,2)*c112 + mab(1,2)*c112t + dmabt(2,2)*c122 + mab(2,2)*c122t);
        
        dq(2) = Eb*(dmab2(1,2) + dmab2(2,2) ...
              + dmabp(1,2)*(c111 + 2*c221) + mab(1,2)*(c111p + 2*c221p) ...
              + dmabp(2,2)*(c112 + 2*c222) + mab(2,2)*(c112p + 2*c222p) ...
              + dmabp(1,1)*c211 + mab(1,1)*c211p + dmabp(2,1)*c221 + mab(2,1)*c221p);

%%      In plane Tension
%       Deformation gradient tensor
        Fd = dxt(:,i,j)*c1R(:,i,j)' + dxp(:,i,j)*c2R(:,i,j)';
        
%       Surface projection operator
        P = eye(3) - nk(:,i,j)*nk(:,i,j)';
        
%       Cauchy-Green
        V2 = Fd*Fd';
        
%       Principal strains
        [ev,lams] = eigs(V2);
        es = sqrt(diag(lams));
%       Normalized eigenvectors        
        ev1 = ev(:,1)/norm(ev(:,1));
        ev2 = ev(:,2)/norm(ev(:,2));
        
%       Strain invariants
        I1 = es(1)^2 + es(2)^2 - 2;
        I2 = es(1)^2*es(2)^2 - 1;
        
%       Covariant in plane deformation gradient tensor !!!! just a check
        Fdv = [dxt(:,i,j)'*Fd*dxt(:,i,j),dxt(:,i,j)'*Fd*dxp(:,i,j)
               dxp(:,i,j)'*Fd*dxt(:,i,j),dxp(:,i,j)'*Fd*dxp(:,i,j)];
%       Contravariant
        Fdn = [c1'*Fd*c1,c1'*Fd*c2
               c2'*Fd*c1,c2'*Fd*c2];
           
%       In plane strain tensor
        eps = (Fdn*Fdv' - eye(2))/2;
        
%       In plane tension
        tau = 0.5*(B/(es(1)*es(2))*(I1 + 1)*V2 ... % Numerical instability here: C increases I2 when it's at ~1e-16 to ~1e-14
            +  es(1)*es(2)*(C*I2-B)*P);
        
%       Principal tensions (not really needed, but could be good to check)
%         taup1 = es(1)/es(2)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);
%         taup2 = es(2)/es(1)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);

%       Matrix in surface coordinates (contravariant)
        tau11 = c1'*tau*c1;
        tau12 = c1'*tau*c2;
        tau21 = c2'*tau*c1;
        tau22 = c2'*tau*c2;
        
%       Now for the derivatives: first derivative of F/V via chain rule
        dFdt = dxt2(:,i,j)*c1R(:,i,j)'  + dxtp(:,i,j)*c2R(:,i,j)' ...
            + dxt(:,i,j)*c1tR(:,i,j)' + dxp(:,i,j)*c2tR(:,i,j)';
        
        dFdp = dxtp(:,i,j)*c1R(:,i,j)'  + dxp2(:,i,j)*c2R(:,i,j)' ...
            + dxt(:,i,j)*c1pR(:,i,j)' + dxp(:,i,j)*c2pR(:,i,j)';
        
        dV2t = dFdt*Fd' + Fd*dFdt';
        dV2p = dFdp*Fd' + Fd*dFdp';
        
%       We need the derivative of the eigenvalues wrt V2 now to get dlamdt
        dlam1dV2 = ev1*ev1';
        dlam2dV2 = ev2*ev2';
        
        des1t = 0.5/sqrt(lams(1,1))*sum(dot(dlam1dV2,dV2t));
        des1p = 0.5/sqrt(lams(1,1))*sum(dot(dlam1dV2,dV2p));
        des2t = 0.5/sqrt(lams(2,2))*sum(dot(dlam2dV2,dV2t));
        des2p = 0.5/sqrt(lams(2,2))*sum(dot(dlam2dV2,dV2p));
        
%       Derivatives of strain invariants
        dI1t = 2*es(1)*des1t + 2*es(2)*des2t;
        dI1p = 2*es(1)*des1p + 2*es(2)*des2p;
        dI2t = (2*es(1)*es(2)^2)*des1t + (es(1)^2*2*es(2))*des2t;
        dI2p = (2*es(1)*es(2)^2)*des1p + (es(1)^2*2*es(2))*des2p;
        
%       Derivatives of projection tensor (double neg in 2nd part)
        dPt = -dnt*-nk(:,i,j)' + nk(:,i,j)*dnt';
        dPp = -dnp*-nk(:,i,j)' + nk(:,i,j)*dnp';
        
%       And finally the big one: derivatives of tau, separated via chain
%       rule, starting with es(1), then es(2)
        dtaut = (-B/(2*es(1)^2*es(2))*(I1+1)*V2 + es(2)/2*(C*I2-B)*P)*des1t ...
              + (-B/(2*es(1)*es(2)^2)*(I1+1)*V2 + es(1)/2*(C*I2-B)*P)*des2t ...
              + 0.5*B/(es(1)*es(2))*V2*dI1t ... % Invariants
              + 0.5*es(1)*es(2)*C*P*dI2t ...
              + 0.5*B/(es(1)*es(2))*(I1+1)*dV2t ...% Tensors
              + 0.5*es(1)*es(2)*(I2-B)*dPt;
          
        dtaup = (-B/(2*es(1)^2*es(2))*(I1+1)*V2 + es(2)/2*(C*I2-B)*P)*des1p ...
              + (-B/(2*es(1)*es(2)^2)*(I1+1)*V2 + es(1)/2*(C*I2-B)*P)*des2p ...
              + 0.5*B/(es(1)*es(2))*V2*dI1p ... % Invariants
              + 0.5*es(1)*es(2)*C*P*dI2p ...
              + 0.5*B/(es(1)*es(2))*(I1+1)*dV2p ...% Tensors
              + 0.5*es(1)*es(2)*(I2-B)*dPp;
          
%       Now put into dtau matrix
        dtauab = [c1'*dtaut*c1, c1'*dtaut*c2;
                  c2'*dtaup*c1, c2'*dtaup*c2];
        
%       Covariant curvature tensor
        bv = [L,M;M,N];
%       Mixed curvature tensor (raise the column)
        bm(1,1) = bv(1,1)*gn(1,1) + bv(1,2)*gn(2,1);
        bm(1,2) = bv(1,1)*gn(1,2) + bv(1,2)*gn(2,2);
        bm(2,1) = bv(2,1)*gn(1,1) + bv(2,2)*gn(2,1);
        bm(2,2) = bv(2,1)*gn(1,2) + bv(2,2)*gn(2,2);
        
%       Covariant divergence of tau in theta, then phi
        cvt = dtauab(1,1) + dtauab(2,1) + tau11*(2*c111 + c221) ...
            + tau21*(2*c112 + c222) + tau12*c112 + tau22*c122;
        cvp = dtauab(1,2) + dtauab(2,2) + tau12*(c111 + 2*c221) ...
            + tau22*(c112 + 2*c222) + tau11*c211 + tau21*c221;
%       Covariant divergence of Q        
        cq = dq(1) + dq(2) + c111*q1 + c112*q2 + c221*q1 + c222*q2;
        
%       And finally the equilibrium equations to get the forces
        fab(1,i,j) = bm(1,1)*q1 + bm(2,1)*q2 - cvt;
        fab(2,i,j) = bm(1,2)*q1 + bm(2,2)*q2 - cvp;
        fab(3,i,j) = -cq - tau11*bv(1,1) - tau12*bv(1,2) - tau21*bv(2,1) - tau22*bv(2,2);
        
%       These are in terms of (contravariant) surface vectors, put them into Cartesian
        myf(:,i,j) = fab(1,i,j)*dxt(:,i,j) + fab(2,i,j)*dxp(:,i,j) + fab(3,i,j)*-nk(:,i,j);
        
        
%       Just to check if derivs are ok
        ggg(i,j) = es(1);
        dggg(i,j) = I2;
        ddgg(i,j) = es(2);
        dddg(i,j) = I1;
    end
end

% Some imaginary components may slip in, so let's make them real
myf = real(myf);

% Harmmonic transforms of the forces. Second spot that introduces some
% inaccuracy, we need this for efficiency purposes, so we don't have to
% calculate the exact force at every point.

% myfS = SpT(Ytf,myf,thf,phf);
% 
% myfS1 = myfS(1,:).';
% myfS2 = myfS(2,:).';
% myfS3 = myfS(3,:).';

myfS1= SpT(Ytf,squeeze(myf(1,:,:)),thf,phf);
myfS2= SpT(Ytf,squeeze(myf(2,:,:)),thf,phf);
myfS3= SpT(Ytf,squeeze(myf(3,:,:)),thf,phf);

% Check energy spectrum of local forces
% fabS1 = SpT(Ytf,squeeze(fab(1,:,:)),thf,phf);
% fabS2 = SpT(Ytf,squeeze(fab(2,:,:)),thf,phf);
% fabS3 = SpT(Ytf,squeeze(fab(3,:,:)),thf,phf);

% Find where the force is distributed on the spectrum to see if aliasing is
% needed really.

for n=0:q
    Est(n+1) = norm(norm([myfS1(n^2+1:(n+1)^2),myfS2(n^2+1:(n+1)^2),myfS3(n^2+1:(n+1)^2)]));
    Es1(n+1) = norm(myfS1(n^2+1:(n+1)^2));
    Es2(n+1) = norm(myfS2(n^2+1:(n+1)^2));
    Es3(n+1) = norm(myfS3(n^2+1:(n+1)^2));
end

%% Calculation of velocity constants via fluid problem

% Calculate the integral
ip = 0;
ic = 0;
A(:) = 0;
b(:) = 0;

% To get just the motion from the stress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bb = b;

% First loop: Inner integrals at Gauss points
for i = 1:nt
    for j = 1:np
        
%       Total Galerkin mode count
        ic = ic+1;
        
%       Velocity at colloc point
        Uc = U + dU*x(:,i,j);
        
%       Rotation Matrix
        t1 = [cos(phi(j)),-sin(phi(j)),0;sin(phi(j)),cos(phi(j)),0;0,0,1];
        t2 = [cos(-tht(i)),0,sin(-tht(i));0,1,0;-sin(-tht(i)),0,cos(-tht(i))];
        t3 = [cos(-phi(j)),-sin(-phi(j)),0;sin(-phi(j)),cos(-phi(j)),0;0,0,1]; 
        Tx = t1*t2*t3;
        
%       Partials of the rotation matrix
        Txt = [-cos(phi(j))^2*sin(tht(i)), -sin(phi(j))*cos(phi(j))*sin(tht(i)), -cos(tht(i))*cos(phi(j));
               -sin(phi(j))*sin(tht(i))*cos(phi(j)), -sin(phi(j))^2*sin(tht(i)), -cos(tht(i))*sin(phi(j));
               cos(tht(i))*cos(phi(j)), cos(tht(i))*sin(phi(j)), -sin(tht(i))];
           
        Txp =[-sin(2*phi(j))*cos(tht(i)) + sin(2*phi(j)), cos(2*phi(j))*cos(tht(i)) - cos(2*phi(j)),  sin(tht(i))*sin(phi(j));
               cos(2*phi(j))*cos(tht(i)) - cos(2*phi(j)), sin(2*phi(j))*cos(tht(i)) - sin(2*phi(j)), -sin(tht(i))*cos(phi(j));
              -sin(tht(i))*sin(phi(j)), sin(tht(i))*cos(phi(j)), 0];
          
%       Rotated harmonic coefficients of geometry
%       !!!!!!!!!! Additionally, I should make it so that SpHRot does the
%       rotation for the three spherical harmonics all at the same time,
%       instead of each individually. Would REALLY cut cost
        x1mnr = SpHRot(x1mn,phi(j),-tht(i),0-phi(j), myfacs);
        x2mnr = SpHRot(x2mn,phi(j),-tht(i),0-phi(j), myfacs);
        x3mnr = SpHRot(x3mn,phi(j),-tht(i),0-phi(j), myfacs);
        
%       Nonrotated locations of rotated grid points.        
        xcg(1,:,:) = SpHReconst(x1mnr,Yt);
        xcg(2,:,:) = SpHReconst(x2mnr,Yt);
        xcg(3,:,:) = SpHReconst(x3mnr,Yt);
        xcg = real(xcg);
        
%       The location of the new pole, should just be location of
%       point at current theta/phi
        xcr = x(:,i,j);
          
%       Gauss points after rotation, in non-rotated reference. 
%       I think the instability could possible be coming from here. We're
%       doing some nonlinear operations on the surface
        for i2 = 1:nt
            for j2 = 1:np
%               Gauss point in rotated reference frame
                gp = [sin(tht(i2))*cos(phi(j2)); sin(tht(i2))*sin(phi(j2)); cos(tht(i2))];
%               We want the rotated coordinates in a nonrotated frame.
%               Tx rotates the colloc point to (0,0,1), so we need the
%               inverse to rotate the GP's. Tx is orthogonal, so just take
%               transpose.
                gp = Tx'*gp;
                phit(i2,j2) = atan2(gp(2),gp(1));
%               Because it can sometimtes be gp(3) - 1 = 2e-16
                if gp(3)> 1; gp(3) =  1;end
                if gp(3)<-1; gp(3) = -1;end
                thet(i2,j2) = acos(gp(3));

%               Using three Spherical harmonics for the three vector coordinates
%               means that the reconstruction will give e.g. the x value in the
%               NON-ROTATED referece frame, USING THE ROTATED reference frame
%               angles (for example, say we rotate the reference 90 deg so the pole
%               is at the x axis in the non rotated frame. Then, an angle of 0,0 in
%               the rotated frame for the x coord returns the value of X at 90,0).

%               Some quick partial derivatives at locations in nonrotated
                dxtg(:,i2,j2) = 0;
                dxpg(:,i2,j2) = 0;
                
                it = 0;
                ih = 0;
                for n = 0:p
                    ih = ih + 1;
%                   Careful here: see comment in ThetDer
                    Yg = Yt{ih};
                    Ypcur = 1i*Yg(:,i2,j2); % To speed up the derivatives
                    im = 0;
                    for m = -n:n
                        im = im+1;
                        it = it+1;
                        
                        td1 = Ytd1_c(i2,j2,it);
                        dxtg(1,i2,j2) = dxtg(1,i2,j2) + x1mnr(it)*td1;
                        dxtg(2,i2,j2) = dxtg(2,i2,j2) + x2mnr(it)*td1;
                        dxtg(3,i2,j2) = dxtg(3,i2,j2) + x3mnr(it)*td1;
                        
                        Ymn = m*Ypcur(im);
                        dxpg(1,i2,j2) = dxpg(1,i2,j2) + x1mnr(it)*Ymn;
                        dxpg(2,i2,j2) = dxpg(2,i2,j2) + x2mnr(it)*Ymn;
                        dxpg(3,i2,j2) = dxpg(3,i2,j2) + x3mnr(it)*Ymn;
                    end
                end
%               Now we can get the Jacobian here, and the normal vector
                
%               Inward normal
                nkg(:,i2,j2) = real(cross(dxtg(:,i2,j2),dxpg(:,i2,j2)));
                nkg(:,i2,j2) = -nkg(:,i2,j2)./norm(nkg(:,i2,j2));
                
%               Jacobian (area element) via fundamental forms
                Jg(i2,j2) = real(sqrt(dot(dxtg(:,i2,j2),dxtg(:,i2,j2)) ... 
                          * dot(dxpg(:,i2,j2),dxpg(:,i2,j2)) ...
                          - dot(dxtg(:,i2,j2),dxpg(:,i2,j2))^2));
            end
        end
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        Ypt = SpHarmT(p,thet,phit,myfacs);
%         YptB = SpHarmT(q,thet,phit,myfacs);
        
%       Forces on rotated grid. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SHOULD
%       JUST BE ABLE TO RECONSTRUCT THESE IN THE LOOP LIKE THE A INTEGRAL
%       This was mewt with massive failure when I tried it though.
        Nmyf(1,:,:) = SpHReconst(myfS1,Ypt);
        Nmyf(2,:,:) = SpHReconst(myfS2,Ypt);
        Nmyf(3,:,:) = SpHReconst(myfS3,Ypt);
        
%       Loop over harmonics
        for n = 0:0 % SPEED MODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%           All spherical harmonics of order n evaluated at (rotated) integ points
            Ypcur = Ypt{n+1};%SpHarm(n,thet,phit);
%           All harmonics at nonrotated points
            Ypcur2 = Yt{n+1};
            
            im = 0;
            for m = -n:n
                ih = ih+1;
                im = im+1;
%               SpHarms order n, degree m eval'd at integ points
                Y = squeeze(Ypcur(im,:,:));
                At(:) = 0;
                for ig = 1:nt
                    for jg = 1:np
                        r = (xcr-xcg(:,ig,jg));
%                         v = Tij(r,nkg(:,ig,jg));
% %                       Integral addition                        
%                         At = At + v*Y(ig,jg)*Jg(ig,jg)*ws(ig)*dphi; % RHS
%                         here can be vectorized
%                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
%                       Only need to calc B once per colloc point
                        if(n==0)
                            r = xcr-xcg(:,ig,jg);
                            v = Gij(r);
                            bt = bt + v*Nmyf(:,ig,jg)*Jg(ig,jg)*ws(ig);
                        end
                    end
                end
%               Plain spherical harmonics addition representing non-integ
%               part (need spherical harmonics at nonrotated Gauss Points)
                Yg = squeeze(Ypcur2(im,:,:));
                
                At = At*(1-lam)/(1+lam) - Yg(i,j)*4*pi*eye(3);
 
                col = 3*(ih-1)+1;
                
%               Integral calc'd for colloc/harm combo. Put in A
%                 A(row:row+2,col:col+2) = A(row:row+2,col:col+2) + At;
%                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end
        end
%       Forcing integral calc'd for a colloc point! Put in B
        bb(row:row+2) = bb(row:row+2) + bt*dphi/mu/(1+lam); %!!!!!!!!!!!!!!!!!!!!!!!!!
        b(row:row+2) = b(row:row+2) + bt*dphi/mu/(1+lam) - Uc*8*pi/(1+lam);
    end
end

% Outer, Galerkin integral (done over a sphere!)
A2(:) = 0;
b2(:) = 0;
it = 0;

bb2 = b2; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Loop over outer product harmonics (a constant value here is a row)
for n = 0:p
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        it = it+1;
        row = 3*it - 2;
        bt(:) = 0;
        bt2 = bt; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Y = squeeze(Ypcur(im,:,:));
        
%       Loop over inner harmonics (a constant value here is a column)
        im2 = 0;
        for n2 = 0:0 % SPEED MODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for m2 = -n2:n2
                im2 = im2+1;
                col = 3*im2 - 2;
%               Loop over Gauss points (these sum to get one row/col value)
                At(:) = 0;
                ii = 0;
                for i = 1:nt
                    for j = 1:np
                        ii = ii+1;
%                       Get value of inner integral at Gauss point
%                         v = A(3*ii-2:3*ii,3*im2-2:3*im2);
%                       Integrate with that value
%                         At = At + v*conj(Y(i,j))*wg(i)*dphi;
                        
%                       Integrate b now (equivalent to finding SpH coeffs
%                       for RHS).
                        if n2 == 0
                            bt = bt + b(3*ii-2:3*ii)*conj(Y(i,j))*wg(i)*dphi;
                            bt2 = bt2 + bb(3*ii-2:3*ii)*conj(Y(i,j))*wg(i)*dphi; %!!!!!!!!!!!
                        end
                    end
                end
                
                A2(row:row+2,col:col+2) = At;
            end
        end
        b2(row:row+2) = bt;
        bb2(row:row+2) = bt2; %!!!!!!!!!!!
    end 
end
% Throw out highest order coefficients (give lots of error, unstable)
% consistent with Zhao, 10.2, paragraph 2, before eq 44, and Rahi, p 777,
% first sentence

b2(3*p^2+1:3*(p+1)^2) = 0;
A2 = eye(3*(p+1)^2)*-4*pi; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ut = A2\b2;
ut(abs(ut)<1e-12) = 0;
%% Convection of points and new surface constants

% Throw out highest order coefficients (give lots of error, unstable)
% consistent with Zhao, 10.2, paragraph 2, before eq 44, and Rahi, p 777,
% first sentence
ut(3*p^2+1:end) = 0;

u1 = ut(1:3:end);
u2 = ut(2:3:end);
u3 = ut(3:3:end);
x1mn = x1mn + u1*dt;
x2mn = x2mn + u2*dt;
x3mn = x3mn + u3*dt;

% Energy Spectra
for n=0:p
    Ex(n+1) = norm(norm([x1mn(n^2+1:(n+1)^2),x2mn(n^2+1:(n+1)^2),x3mn(n^2+1:(n+1)^2)]));
    Eu(n+1) = norm(norm([u1(n^2+1:(n+1)^2),u2(n^2+1:(n+1)^2),u3(n^2+1:(n+1)^2)]));
end

% Save the shapes for post-processing
xmns(:,1,cts+1) = x1mn;
xmns(:,2,cts+1) = x2mn;
xmns(:,3,cts+1) = x3mn;

% Cartesian displacements
ua11(1,:,:) = real(SpHReconst(u1,Yt));
ua11(2,:,:) = real(SpHReconst(u2,Yt));
ua11(3,:,:) = real(SpHReconst(u3,Yt));


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! switching vel here to just stress
ut = A2\bb2;

u1 = ut(1:3:end);
u2 = ut(2:3:end);
u3 = ut(3:3:end);
ua(1,:,:) = real(SpHReconst(u1,Yt));
ua(2,:,:) = real(SpHReconst(u2,Yt));
ua(3,:,:) = real(SpHReconst(u3,Yt));

% Plots
clf;
x1 = real(SpHReconst(x1mn,Yr));
x2 = real(SpHReconst(x2mn,Yr));
x3 = real(SpHReconst(x3mn,Yr));

% New x's
x(1,:,:) = real(SpHReconst(x1mn,Yt));
x(2,:,:) = real(SpHReconst(x2mn,Yt));
x(3,:,:) = real(SpHReconst(x3mn,Yt));

fpl = zeros(nt,np);
for i = 1:nt
    for j =1:np
        fpl(i,j) = norm(myf(1:2,i,j));
    end
end

% Force reconstructys
myf21 = real(SpHReconst(myfS1,Yr));
myf22 = real(SpHReconst(myfS2,Yr));
myf23 = real(SpHReconst(myfS3,Yr));
% tau11r = real(SpHReconst(f11,Yr));
% tau22r = real(SpHReconst(f22,Yr));
% tau21r = real(SpHReconst(f21,Yr));
myft = sqrt(myf21.^2 + myf22.^2 + myf23.^2);

% Integral of f over surface
fintx = real(bb(1:3:end));
fintx = reshape(fintx,np,nt)';
finty = real(bb(2:3:end));
finty = reshape(finty,np,nt)';
fintz = real(bb(3:3:end));
fintz = reshape(fintz,np,nt)';

% Quantities to keep track of
trcpnt(cts) = real(SpHReconst(x1mn,Ytrc));%x(1,floor(nt/2)+1,1);
trcvel(cts) = real(SpHReconst(u1,Ytrc));%ua(1,floor(nt/2)+1,1);
trcpntz(cts) = real(SpHReconst(x3mn,Ytrc));%x(1,floor(nt/2)+1,1);
trcvelz(cts) = real(SpHReconst(u3,Ytrc));%ua(3,floor(nt/2)+1,1);

% q=quiver3(reshape(x1,[1,numel(x1)]),reshape(x2,[1,numel(x1)]),reshape(x3,[1,numel(x1)]),reshape(myf21,[1,numel(x1)]),reshape(myf22,[1,numel(x1)]),reshape(myf23,[1,numel(x1)]));

h1 = subplot(2,1,1);
set(h1, 'Units', 'normalized');
set(h1, 'Position', [-.1, 0.5, 1.15, .6]);

surf(x1,x2,x3,myf23,'edgecolor','none','FaceColor',[1 0 0], ...
     'FaceAlpha',0.75,'FaceLighting','gouraud')
lightangle(gca,150,50)
set(gca,'nextplot','replacechildren','visible','off')
% Top down
% view(0,90);
% Side
view(0,0);
axis([-2,2,0,2,-2,2])
pbaspect([1,.5,1])
hold on
scatter3(trcpnt(cts),real(SpHReconst(x2mn,Ytrc)),real(SpHReconst(x3mn,Ytrc)),75,'go','filled');

% Just from stress - magenta
% quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua(1,:,:),[1,numel(x(1,:,:))]),reshape(ua(2,:,:),[1,numel(x(1,:,:))]),reshape(ua(3,:,:),[1,numel(x(1,:,:))]),'m')
% Force - green
% quiver3(reshape(xf(1,:,:),[1,numel(xf(1,:,:))]),reshape(xf(2,:,:),[1,numel(xf(1,:,:))]),reshape(xf(3,:,:),[1,numel(xf(1,:,:))]),reshape(real(myf(1,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(2,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(3,:,:)),[1,numel(xf(1,:,:))]),'g')
% Stress and fluid - blue
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua11(1,:,:),[1,numel(x(1,:,:))]),reshape(ua11(2,:,:),[1,numel(x(1,:,:))]),reshape(ua11(3,:,:),[1,numel(x(1,:,:))]),'b')

h2 = subplot(2,1,2);
set(h2, 'Units', 'normalized');
set(h2, 'Position', [0.05, 0, 1, .7]);

surf(x1,x2,x3,myf23,'edgecolor','none','FaceColor',[1 0 0], ... %!!!!!!!!!!!!! make full cell?
      'FaceAlpha',0.75,'FaceLighting','gouraud')
lightangle(gca,150,50)
set(gca,'nextplot','replacechildren','visible','off')
set(gcf, 'color', 'white');
% Top down
% view(0,90);
% Side
view(45,45);
axis([-2,2,0,2,-2,2])
pbaspect([1,.5,1])
hold on
scatter3(trcpnt(cts),real(SpHReconst(x2mn,Ytrc)),real(SpHReconst(x3mn,Ytrc)),75,'go','filled');


% Just from stress - magenta, includes highest order!
% quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua(1,:,:),[1,numel(x(1,:,:))]),reshape(ua(2,:,:),[1,numel(x(1,:,:))]),reshape(ua(3,:,:),[1,numel(x(1,:,:))]),'m')
% Force - green
% quiver3(reshape(xf(1,:,:),[1,numel(xf(1,:,:))]),reshape(xf(2,:,:),[1,numel(xf(1,:,:))]),reshape(xf(3,:,:),[1,numel(xf(1,:,:))]),reshape(real(myf(1,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(2,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(3,:,:)),[1,numel(xf(1,:,:))]),'g')
% Stress and fluid - blue
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua11(1,:,:),[1,numel(x(1,:,:))]),reshape(ua11(2,:,:),[1,numel(x(1,:,:))]),reshape(ua11(3,:,:),[1,numel(x(1,:,:))]),'b')


% %vel plots
% res = 10;
% drs = 4/res;
% xqv1 = zeros(1,res*res);
% xqv2 = xqv1;
% xqv3 = xqv1;
% vqv1 = xqv1;
% vqv2 = xqv1;
% vqv3 = xqv1;
% ntt=0;
% for i = 1:res
%     for j = 1:res
%         ntt = ntt+1;
%         xqv1(ntt) = -2 + (i-1)*drs;
%         xqv3(ntt) = -2 + (j-1)*drs;
%         vt = U + dU*[xqv1(ntt);xqv2(ntt);xqv3(ntt)];
%         vqv1(ntt) = vt(1);
%         vqv2(ntt) = vt(2);
%         vqv3(ntt) = vt(3);
%     end
% end
% quiver3(xqv1,xqv2,xqv3,vqv1,vqv2,vqv3);


% gif maker
drawnow
% Capture the plot as an image 
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(frame.cdata,256,'nodither');
% Write to the GIF File 
if cts == 1 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end 

% some output
disp([max(abs(x3mn-ix3mn)), real(u3(1)), max(max(squeeze(ua(1,:,:)))),max(max(max(abs(myf(:,:,:))))),cts]);

 end

















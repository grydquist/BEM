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
B = 0.005;
C = 100;
% B=1;
% C=1;

% Bending modulus
Eb = 1;%0.0669; % Non-dimensional

% Total time steps
NT = 200;
% dt = 0.00001;
dt = 0.0001;

% Velocity and gradient
U = [0;0;0];
dU = [0,5,0;5,0,0;0,0,0];
% dU(:) = 0;
%% Reference Shape Definition

% Order of the force and velocity transforms
p = 3;

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

% Weights for rotated integral calculation
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:gf*p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmT(p,th,ph);

% Spherical harmonic coefficients defining the shape of ref & def:
% xmnR = [2*sqrt(pi),0,0,0,0,0,0,0,0];%,0,0,0,0,0,0,0,0];
% xmnR = [2*sqrt(pi),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% xmnR = [2*sqrt(pi),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
xmnR = zeros((p+1)^2,1);
xmnR(1) = 2*sqrt(pi);
% xmnR = [4.5464,0,0,0,0,.1i,-0.6582,.1i,0];

% Get harmonics for the individual parts of the position vector (would be
% good to do this for the reference too, but don't need to at this point)
rrr = SpHReconst(xmnR,Yt);
x1 = rrr.*sin(th).*cos(ph);
x2 = rrr.*sin(th).*sin(ph);
x3 = rrr.*cos(th);

x1mn = SpT(Yt,x1,th,ph);
x2mn = SpT(Yt,x2,th,ph);
x3mn = SpT(Yt,x3,th,ph);

% 
% x1mn = zeros((p+1)^2,1);
% x3mn = x1mn;
% 
% x1mn(1:4) = [0,2*sqrt(2*pi/3),0,-2*sqrt(2*pi/3)]/2;
% x2mn = 1i*abs(x1mn);
% x3mn(1:4) = [0,0,4*sqrt(pi/3),0]/2;


% Evaluation of the needed geometrical infor for the reference state
[xR, nR, JR, dxtR, dxpR, kR, kdR] = SpHDer(xmnR, Yt, th,ph);
% Contravariant form of reference tangent

c1R = zeros(3,nt,np);
c2R = c1R;
for i = 1:nt
    for j = 1:np
        c1R(:,i,j) = cross(dxpR(:,i,j),-nR(:,i,j))/JR(i,j);
        c2R(:,i,j) = cross(-nR(:,i,j),dxtR(:,i,j))/JR(i,j);
    end
end

%% Preallocation of all the solid things I need

% Actual locations
x = zeros(3,nt,np);
rr = zeros(nt,np);

% Total derivatives, in Cartesian
dxp = zeros(3,nt,np);
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

% Normal Vector (inward)
nk = zeros(3,nt,np);

% Fundamental forms/metric tensor
E = zeros(nt,np);
F = E;
G = E;
J = E;
L = E;
M = E;
N = E;
k = E;

% Containers for the stresses
tau11 = zeros(nt,np);
tau12 = tau11;
tau21 = tau11;
tau22 = tau11;
q1 = tau11;
q2 = tau11;
dmab = zeros(2,2);

% Saving the Christoffels
c111 = zeros(nt,np);
c112 = c111;
c122 = c111;
c222 = c111;
c221 = c111;
c211 = c111;

dtau  = zeros(2,2);
dq = zeros(2,1);
bm = dtau;
bn = bm;

myf = zeros(3,nt,np);
Nmyf = myf;
fab = zeros(3,1);
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
dxtg = zeros(3,1);
dxpg = dxtg;
xcg = zeros(3,nt,np);
ua = zeros(3,nt,np);
% Harmonics at north pole
Ytcr= SpHarmT(p,0,0);

%% Time stepping loop
for cts = 1:NT
%% Calculation of residual force at convected points
% After calculation, do SpHarm transform
if(cts==2);dU(:)=0;end %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Get Geometric information for current configuration (Perhaps the best way
% to do this is with a subroutine that takes in the SpHarm consts for the
% def, and infor for the ref, and material constants and just spits out the
% forces. In FORTRAN this information would be packaged up in objects, and
% we could just pass pointers to the objects!)

ih = 0;
it = 0;
x(:) = 0;
dxt(:) = 0;
dxp(:) = 0;
dxt2(:) = 0;
dxp2(:) = 0;
dxtp(:) = 0;
dxt3(:) = 0;
dxp3(:) = 0;
dxt2p(:) = 0;
dxtp2(:) = 0;

% Calcaulating surface derivatives
% Loop over harmonics
for n = 0:p
    ih = ih+1;
%   This is assuming that the deformed and reference are at the same tht
%   and phi, which they should be
    Ypcur = Yt{ih};
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
        td1 = ThetDer(Yt,ph,n,m,1);
        
%       Distance from zero
        x(1,:,:) = squeeze(x(1,:,:)) + f1*Ym;
        x(2,:,:) = squeeze(x(2,:,:)) + f2*Ym;
        x(3,:,:) = squeeze(x(3,:,:)) + f3*Ym;
        x = real(x);
        
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
        td2 = ThetDer(Yt,ph,n,m,2);

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
        td3 = ThetDer(Yt,ph,n,m,3);
        dxt3(1,:,:) = squeeze(dxt3(1,:,:)) + f1*td3;
        dxt3(2,:,:) = squeeze(dxt3(2,:,:)) + f2*td3;
        dxt3(3,:,:) = squeeze(dxt3(3,:,:)) + f3*td3;
    end
end


for i = 1:nt
    for j = 1:np
        
%       Normal vector (inward)
        nk(:,i,j) = cross(dxt(:,i,j),dxp(:,i,j));
        nk(:,i,j) = -nk(:,i,j)./norm(nk(:,i,j));
        nk = real(nk);
        
%       Radial displacement
        rr(i,j) = sqrt(x(1,i,j).^2 + x(2,i,j).^2 + x(3,i,j).^2);
        
%       Jacobian via fundamental forms
        E(i,j) = dot(dxt(:,i,j),dxt(:,i,j));
        F(i,j) = dot(dxt(:,i,j),dxp(:,i,j));
        G(i,j) = dot(dxp(:,i,j),dxp(:,i,j));
        J(i,j) = sqrt(E(i,j).*G(i,j) - F(i,j).^2);

%       Calculate mean curvature with some fundamental forms (appendices in
%       Veeranpani paper)
        L(i,j) = dot(dxt2(:,i,j),-nk(:,i,j));
        M(i,j) = dot(dxtp(:,i,j),-nk(:,i,j));
        N(i,j) = dot(dxp2(:,i,j),-nk(:,i,j));
        k(i,j) = (E(i,j)*N(i,j) - 2*F(i,j)*M(i,j) + G(i,j)*L(i,j))/(2*J(i,j)^2);
        
%       First let's do all the cross products we need
        ct2_p = cross(dxt2(:,i,j),dxp(:,i,j));
        ct_tp = cross(dxt (:,i,j),dxtp(:,i,j));
        ctp_p = cross(dxtp(:,i,j),dxp(:,i,j));
        ct_p2 = cross(dxt (:,i,j),dxp2(:,i,j));
        
%       Jacobian Partials        
        Jt = dot(-nk(:,i,j),ct2_p + ct_tp);
        Jp = dot(-nk(:,i,j),ctp_p + ct_p2);

%       Normal vector partials
        dnt = (1/J(i,j))*(ct2_p + ct_tp - Jt*-nk(:,i,j));
        dnp = (1/J(i,j))*(ctp_p + ct_p2 - Jp*-nk(:,i,j));
     
%       Fundamental form partials        
        Et = 2*dot(dxt2(:,i,j),dxt(:,i,j));
        Ep = 2*dot(dxtp(:,i,j),dxt(:,i,j));
        Gt = 2*dot(dxtp(:,i,j),dxp(:,i,j));
        Gp = 2*dot(dxp (:,i,j),dxp2(:,i,j));
        Ft = dot(dxt2(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxtp(:,i,j));
        Fp = dot(dxtp(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxp2(:,i,j));
        Lt = dot(dxt3(:,i,j),-nk(:,i,j))  + dot(dxt2(:,i,j),dnt);
        Lp = dot(dxt2p(:,i,j),-nk(:,i,j)) + dot(dxt2(:,i,j),dnp);
        Nt = dot(dxtp2(:,i,j),-nk(:,i,j)) + dot(dxp2(:,i,j),dnt);
        Np = dot(dxp3(:,i,j),-nk(:,i,j))  + dot(dxp2(:,i,j),dnp);
        Mt = dot(dxt2p(:,i,j),-nk(:,i,j)) + dot(dxtp(:,i,j),dnt);
        Mp = dot(dxtp2(:,i,j),-nk(:,i,j)) + dot(dxtp(:,i,j),dnp);

%       Curvature derivatives        
        kt = -2*Jt./J(i,j).*k(i,j) + 0.5*J(i,j).^-2. ... 
                * (Et.*N(i,j) + Nt.*E(i,j) - 2*Ft.*M(i,j) - 2*Mt.*F(i,j) + Gt.*L(i,j) + Lt.*G(i,j));
        kp = -2*Jp./J(i,j).*k(i,j) + 0.5*J(i,j).^-2. ... 
                * (Ep.*N(i,j) + Np.*E(i,j) - 2*Fp.*M(i,j) - 2*Mp.*F(i,j) + Gp.*L(i,j) + Lp.*G(i,j));       
            
%       Shear transverse
%       Contravariant basis vectors        
        c1 = cross(dxp(:,i,j),-nk(:,i,j))/J(i,j);
        c2 = cross(-nk(:,i,j),dxt(:,i,j))/J(i,j);
%       !!!!!!!!!! HACKY FIRST TIME STEP
        if(cts == 1)
            c1R(:,i,j) = c1;
            c2R(:,i,j) = c2;
            kR(i,j) = k(i,j);
            kdR(1,i,j) = kt;
            kdR(2,i,j) = kp;
        end
        
%       Metric tensor
        g = [E(i,j),F(i,j);F(i,j),G(i,j)];
%       Contravariant g
        gn = inv(g);
%       Partials of COTNRAVARIANT metric tensor
        dgt = -gn*[Et,Ft;Ft,Gt]*gn;
        dgp = -gn*[Ep,Fp;Fp,Gp]*gn;
        
%       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
        c111(i,j) = dot(dxt2(:,i,j),c1);
        c112(i,j) = dot(dxtp(:,i,j),c1);
        c122(i,j) = dot(dxp2(:,i,j),c1);
        c222(i,j) = dot(dxp2(:,i,j),c2);
        c221(i,j) = dot(dxtp(:,i,j),c2);
        c211(i,j) = dot(dxt2(:,i,j),c2);
        
%       Moments and their needed partials
        mab = (k(i,j) - kR(i,j))*gn;
        dmab(1,1) = (kt - kdR(1,i,j))*gn(1,1) + (k(i,j) - kR(i,j))*dgt(1,1);
        dmab(1,2) = (kt - kdR(1,i,j))*gn(1,2) + (k(i,j) - kR(i,j))*dgt(1,2);
        dmab(2,1) = (kp - kdR(2,i,j))*gn(2,1) + (k(i,j) - kR(i,j))*dgp(2,1);
        dmab(2,2) = (kp - kdR(2,i,j))*gn(2,2) + (k(i,j) - kR(i,j))*dgp(2,2);
        
%       Transverse Shears (Covariant divergence of moment tensor)
        q1(i,j) = Eb*(dmab(1,1) + dmab(2,1) + mab(1,1)*(2*c111(i,j) + c221(i,j)) ...
                + mab(2,1)*(2*c112(i,j) + c222(i,j)) + mab(1,2)*c112(i,j) + mab(2,2)*c122(i,j));
        q2(i,j) = Eb*(dmab(1,2) + dmab(2,2) + mab(1,2)*(c111(i,j) + 2*c221(i,j)) ...
                + mab(2,2)*(c112(i,j) + 2*c222(i,j)) + mab(1,1)*c211(i,j) + mab(2,1)*c221(i,j));
            
%%      In plane Tension
%       Deformation gradient tensor
        Fd = dxt(:,i,j)*c1R(:,i,j)' + dxp(:,i,j)*c2R(:,i,j)';
        
%       Surface projection operator
        P = eye(3) - nk(:,i,j)*nk(:,i,j)';
        
%       Surface deformation gradient tensor
        Fs = P*Fd;
        
%       Cauchy-Green
        V2 = Fd*Fd';
        
%       Principal strains
        [ev,lams] = eigs(V2);
        es = sqrt(diag(lams));
        
%       Strain invariants
        I1 = es(1)^2 + es(2)^2 - 2;
        I2 = es(1)^2*es(2)^2 - 1;
        
%       In plane tension
        tau = B/(2*es(1)*es(2))*(I1 + 1)*V2 ...
            +  0.5*es(1)*es(2)*(C*I2-B)*P;
        
%       Principal tensions (not really needed, but could be good to check)
%         taup1 = es(1)/es(2)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);
%         taup2 = es(2)/es(1)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);

%       Matrix in surface coordinates (contravariant)
        tau11(i,j) = c1'*tau*c1;
        tau12(i,j) = c1'*tau*c2;
        tau21(i,j) = c2'*tau*c1;
        tau22(i,j) = c2'*tau*c2;
        
    end
end
% Do the Spherical Harmonic transform so we can take a surface divergence to get the surface stresses

% Spherical harmonic coefficients for each component of the stress
% Pretty inefficient... This is also the first spot where inaccuracy comes
% into the equation. Before, the stresses were calculated exactly for a
% given spherical harmonic representation of the reference and deformed
% state

f11 = SpT(Yt,tau11,th,ph);
f12 = SpT(Yt,tau12,th,ph);
f21 = SpT(Yt,tau21,th,ph);
f22 = SpT(Yt,tau22,th,ph);
fq1 = SpT(Yt,q1,th,ph);
fq2 = SpT(Yt,q2,th,ph);

for i = 1:nt
    for j =1:np
%       Calculate the partials of tau and q        
        dtau(:) = 0;
        dq(:) = 0;
        ih = 0;
        it = 0;
        for k = 0:p
            ih = ih+1; 
            Ypcur = Yt{ih};
            im = 0;
            for l = -k:k
                it = it+1;
                im = im+1;
%               Current values of harmonic at degree and order
                Ym = squeeze(Ypcur(im,:,:));
                
%               Theta partial
                td1 = ThetDer(Yt,ph,k,l,1);
                dtau(1,1) = dtau(1,1) + td1(i,j)*f11(it);
                dtau(1,2) = dtau(1,2) + td1(i,j)*f12(it);
                dq(1) = dq(1) + td1(i,j)*fq1(it);
                
%               Phi partial
                dtau(2,1) = dtau(2,1) + f21(it)*1i*l*Ym(i,j);
                dtau(2,2) = dtau(2,2) + f22(it)*1i*l*Ym(i,j);
                dq(2) = dq(2) + fq2(it)*1i*l*Ym(i,j);
            end
        end
%       Metric tensor
        g = [E(i,j),F(i,j);F(i,j),G(i,j)];
%       Contravariant g
        gn = inv(g);
        
%       Covariant curvature tensor
        bv = [L(i,j),M(i,j);M(i,j),N(i,j)];
%       Mixed curvature tensor (raise the column)
        bm(1,1) = bv(1,1)*gn(1,1) + bv(1,2)*gn(2,1);
        bm(1,2) = bv(1,1)*gn(1,2) + bv(1,2)*gn(2,2);
        bm(2,1) = bv(2,1)*gn(1,1) + bv(2,2)*gn(2,1);
        bm(2,2) = bv(2,1)*gn(1,2) + bv(2,2)*gn(2,2);
        
%       Covariant divergence of tau in theta, then phi
        cvt = dtau(1,1) + dtau(2,1) + tau11(i,j)*(2*c111(i,j) + c221(i,j)) ...
            + tau21(i,j)*(2*c112(i,j) + c222(i,j)) + tau12(i,j)*c112(i,j) + tau22(i,j)*c122(i,j);
        cvp = dtau(1,2) + dtau(2,2) + tau12(i,j)*(c111(i,j) + 2*c221(i,j)) ...
            + tau22(i,j)*(c112(i,j) + 2*c222(i,j)) + tau11(i,j)*c211(i,j) + tau21(i,j)*c221(i,j);
%       Covariant divergence of Q        
        cq = dq(1) + dq(2) + c111(i,j)*q1(i,j) + c112(i,j)*q2(i,j) + c221(i,j)*q1(i,j) + c222(i,j)*q2(i,j);
        
%       And finally the equilibrium equations to get the forces
        fab(1) = bm(1,1)*q1(i,j) + bm(2,1)*q2(i,j) - cvt;
        fab(2) = bm(1,2)*q1(i,j) + bm(2,2)*q2(i,j) - cvp;
        fab(3) = -cq - tau11(i,j)*bv(1,1) - tau12(i,j)*bv(1,2) - tau21(i,j)*bv(2,1) - tau22(i,j)*bv(2,2);
        
%       These are in terms of (contravariant) surface vectors, put them into Cartesian
        myf(:,i,j) = fab(1)*dxt(:,i,j) + fab(2)*dxp(:,i,j) + fab(3)*-nk(:,i,j);
    end
end

% Some imaginary components may slip in, so let's make them real
myf = real(myf);

% Harmmonic transforms of the forces. Second spot that introduces some
% inaccuracy, we need this for efficiency purposes, so we don't have to
% calculate the exact force at every point.
myfS1= SpT(Yt,squeeze(myf(1,:,:)),th,ph);
myfS2= SpT(Yt,squeeze(myf(2,:,:)),th,ph);
myfS3= SpT(Yt,squeeze(myf(3,:,:)),th,ph);

%% Calculation of velocity constants via fluid problem

% Calculate the integral
ip = 0;
ic = 0;
A(:) = 0;
b(:) = 0;

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
        x1mnr = SpHRot(x1mn,phi(j),-tht(i),0-phi(j));
        x2mnr = SpHRot(x2mn,phi(j),-tht(i),0-phi(j));
        x3mnr = SpHRot(x3mn,phi(j),-tht(i),0-phi(j));
        
%       Nonrotated locations of rotated grid points.        
        xcg(1,:,:) = SpHReconst(x1mnr,Yt);
        xcg(2,:,:) = SpHReconst(x2mnr,Yt);
        xcg(3,:,:) = SpHReconst(x3mnr,Yt);
        xcg = real(xcg);
          
%       Gauss points after rotation, in non-rotated reference. !!! likely a better way to do this
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

%       Using three Spherical harmonics for the three vector coordinates
%       means that the reconstruction will give e.g. the x value in the
%       NON-ROTATED referece frame, USING THE ROTATED reference frame
%       angles (for example, say we rotate the reference 90 deg so the pole
%       is at the x axis in the non rotated frame. Then, an angle of 0,0 in
%       the rotated frame for the x coord returns the value of X at 90,0).
%       That's not what we want: we want the x value in the rotated
%       reference frame, when viewed from that reference frame.
%               Some quick partial derivatives at locations in nonrotated
                dxtg(:) = 0;
                dxpg(:) = 0;
                
%               Another inefficieny: we already have these all calculated, but we only want to do ThetDer at one point.                
                Yg = SpHarmT(p,tht(i2),phi(j2));
                it = 0;
                ih = 0;
                for n = 0:p
                    ih = ih + 1;
                    Ypcur = Yg{ih};
                    im = 0;
                    for m = -n:n
                        im = im+1;
                        it = it+1;
                        
                        td1 = ThetDer(Yg,phi(j2),n,m,1);
                        dxtg(1) = dxtg(1) + x1mnr(it)*td1;
                        dxtg(2) = dxtg(2) + x2mnr(it)*td1;
                        dxtg(3) = dxtg(3) + x3mnr(it)*td1;
                        
                        dxpg(1) = dxpg(1) + x1mnr(it)*1i*m*Ypcur(im);
                        dxpg(2) = dxpg(2) + x2mnr(it)*1i*m*Ypcur(im);
                        dxpg(3) = dxpg(3) + x3mnr(it)*1i*m*Ypcur(im);
                    end
                end
%               Now we can get the Jacobian here, and the normal vector
                
%               Getting location of rotated grid points in rotated frame
                xcg(:,i2,j2) = Tx*xcg(:,i2,j2);
                
%               Partials in rotated frame (needed)?
%                 dxtg = Tx*dxtg;
%                 dxpg = Tx*dxpg;

%               Inward normal
                nkg(:,i2,j2) = real(cross(dxtg,dxpg));
                nkg(:,i2,j2) = -nkg(:,i2,j2)./norm(nkg(:,i2,j2));
                nkg(:,i2,j2) = Tx*nkg(:,i2,j2);
%               Jacobian (area element) via fundamental forms
                Jg(i2,j2) = real(sqrt(dot(dxtg,dxtg)*dot(dxpg,dxpg) - dot(dxtg,dxpg)^2));
                
            end
        end
        
%       The location of the new pole (just the radial displacement of the 
%       Gauss point we rotated to)

        xcr = Tx*[SpHReconst(x1mnr,Ytcr),SpHReconst(x2mnr,Ytcr),SpHReconst(x3mnr,Ytcr)]';
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        Ypt = SpHarmT(p,thet,phit);
        Nmyf(1,:,:) = SpHReconst(myfS1,Ypt);
        Nmyf(2,:,:) = SpHReconst(myfS2,Ypt);
        Nmyf(3,:,:) = SpHReconst(myfS3,Ypt);
        
%       Loop over harmonics
        for n = 0:p
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
                        r = Tx'*(xcr-xcg(:,ig,jg));
                        v = Tij(r,Tx'*nkg(:,ig,jg));
%                       Integral addition                        
                        At = At + v*Y(ig,jg)*Jg(ig,jg)*ws(ig)*dphi;
                        
%                       Only need to calc B once per colloc point
                        if(n==0)
                            r = (xcr-xcg(:,ig,jg));
                            v = Gij(r);
%                           Do I also need to rotate the force and make neg???
%                             bt = bt + v*Tx'*-myf(:,ig,jg)*Jg(ig,jg)*ws(ig);
                            bt = bt + v*Nmyf(:,ig,jg)*Jg(ig,jg)*ws(ig);
%                             bt = bt + v*Tx'*Nmyf(:,ig,jg)*Jg(ig,jg)*ws(ig);
                        end
                    end
                end
%               Plain spherical harmonics addition representing non-integ
%               part (need spherical harmonics at nonrotated Gauss Points)
                Yg = squeeze(Ypcur2(im,:,:));
                
                At = At*(1-lam)/(1+lam) - Yg(i,j)*4*pi*eye(3);
 
                col = 3*(ih-1)+1;
                
%               Integral calc'd for colloc/harm combo. Put in A
                A(row:row+2,col:col+2) = A(row:row+2,col:col+2) + At;
            end
        end
%       Forcing integral calc'd for a colloc point! Put in B
        b(row:row+2) = b(row:row+2) + bt*dphi - 4*pi*Uc*2/(1+lam);
    end
end

% Outer, Galerkin integral (done over a sphere!)
A2(:) = 0;
b2(:) = 0;
it = 0;

% Loop over outer product harmonics (a constant value here is a row)
for n = 0:p
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        it = it+1;
        row = 3*it - 2;
        bt(:) = 0;
        Y = squeeze(Ypcur(im,:,:));
        
%       Loop over inner harmonics (a constant value here is a column)
        im2 = 0;
        for n2 = 0:p
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
                        v = A(3*ii-2:3*ii,3*im2-2:3*im2);
%                       Integrate with that value
                        At = At + v*conj(Y(i,j))*wg(i)*dphi;
                        
%                       Integrate b now
                        if n2 == 0
                            bt = bt + b(3*ii-2:3*ii)*conj(Y(i,j))*wg(i)*dphi;
                        end
                    end
                end
                
                A2(row:row+2,col:col+2) = At;
            end
        end
        b2(row:row+2) = bt;
    end
end

% [Us,S,V] = svd(A2);
% Si = 1./S;
% % % Si(Si>1e10)=0;
% % % Si(end) = 0;
% Ai = V*Si*Us';
ut = A2\b2;%Ai*b2;
%% Convection of points and new surface constants
% Including getting new angles for convected points, then go back to the
% top

u1 = ut(1:3:end);
u2 = ut(2:3:end);
u3 = ut(3:3:end);
x1mn = x1mn + u1*dt;
x2mn = x2mn + u2*dt;
x3mn = x3mn + u3*dt;

% Get Cartesian displacements
ind = 0;
ua(:) = 0;
for n = 0:p
%   All spherical harmonics of order n evaluated at integ points
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        Y = Ypcur(im,:,:);
        for d = 1:3
            ind = ind+1;
%           SpHarms order n, degree m eval'd at integ points
            ua(d,:,:) = ua(d,:,:) + Y*ut(ind);
        end
    end
end
ua = real(ua); 

% Displace the material points, and get the new "angles" associated with
% that displacement

% Note that this may not work after time step 1, as ua is calcd at GP's
% x = x + ua*dt;

% Plots
clf;
tmpt = linspace(0,pi,100);tmpp = linspace(0,2*pi,100);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmT(p,tmpth,tmpph);
x1 = real(SpHReconst(x1mn,Yr));
x2 = real(SpHReconst(x2mn,Yr));
x3 = real(SpHReconst(x3mn,Yr));
try1 = sqrt(x1.^2+x2.^2 +x3.^2);
% try1 = squeeze(sqrt(x(1,:,:).^2+x(2,:,:).^2 +x(3,:,:).^2));
% tmpph = squeeze(atan2(x(2,:,:),x(1,:,:)));
% tmpth = squeeze(atan2(sqrt(x(1,:,:).^2 + x(2,:,:).^2),x(3,:,:)));
[Xm,Ym,Zm] = sph2cart(tmpph, tmpth-pi/2, try1);

fpl = zeros(nt,np);
for i = 1:nt
    for j =1:np
        fpl(i,j) = norm(myf(1:2,i,j));
    end
end

% myf21 = SpHReconst(myfS1,Yr);
% myf22 = SpHReconst(myfS2,Yr);
% myf23 = SpHReconst(myfS3,Yr);
% q=quiver3(reshape(x1,[1,numel(x1)]),reshape(x2,[1,numel(x1)]),reshape(x3,[1,numel(x1)]),reshape(myf21,[1,numel(x1)]),reshape(myf22,[1,numel(x1)]),reshape(myf23,[1,numel(x1)]));

surf(Xm,Ym,Zm,'edgecolor','none')
view(0,90);
axis([-2,2,-2,2,-2,2])
pbaspect([1,1,1])
hold on
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua(1,:,:),[1,numel(x(1,:,:))]),reshape(ua(2,:,:),[1,numel(x(1,:,:))]),reshape(ua(3,:,:),[1,numel(x(1,:,:))]))
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(myf(1,:,:),[1,numel(x(1,:,:))]),reshape(myf(2,:,:),[1,numel(x(1,:,:))]),reshape(myf(3,:,:),[1,numel(x(1,:,:))]))
drawnow
disp([real(u1(2)),real(u2(2)), real(u3(1)),cts]);
% disp(x1mn(2))
end






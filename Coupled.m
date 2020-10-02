% Differences: get the deformed state info first. With that, calculate the
% velocity. Only need to get material constants once.
%% Material Constants
% Viscosity outside
mu = 1;
% Viscosity ratio
lam = 1/5;

% Deformation resistance constants
B = 0.005;
C = 100;

% Bending modulus
Eb = 1;%0.0669; % Non-dimensional

% Total time steps
NT = 40;
dt = 0.00001;

% Velocity and gradient
U = [0;0;0];
dU = [0,0,0;0,0,0;0,0,0];
%% Reference Shape Definition

% Order of the force and velocity transforms
p = 1;

% How many coefficients per dimension?
ftot = (p+1)^2;

% Number of points times more than needed for the integrals
gf = 2;

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

% Spherical harmonic coefficients defining the shape of ref & def
xmnR = 2*sqrt(pi);
xmn = zeros(1,ftot);
xmn(1) = 4*sqrt(pi); % In most cases this should just be xmnR


% Order of geometric parameterization
pxmn = sqrt(length(xmnR))-1;
% Check that this is an integer
if floor(pxmn) ~= pxmn
    error('Geometric parameterization has the wrong number of constants')
end

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmT(p,th,ph);

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

% Actual distance
r = zeros(nt,np);
% Actual locations
x = zeros(3,nt,np);

% Derivatives of just harmonics at phi/tht
rt = zeros(nt,np);
rp = rt;

%   Higher order derivatives
% 2nd
rt2 = rt;
rp2 = rt;
rtp = rt;

% 3rd
rt3 = rt;
rp3 = rt;
rt2p= rt;
rtp2= rt;

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
ua = zeros(3,nt,np);

%% Time stepping loop
for ct = 1:NT
%% Calculation of residual force at convected points
% After calculation, do SpHarm transform


% Get Geometric information for current configuration (Perhaps the best way
% to do this is with a subroutine that takes in the SpHarm consts for the
% def, and infor for the ref, and material constants and just spits out the
% forces. In FORTRAN this information would be packaged up in objects, and
% we could just pass pointers to the objects!)

ih = 0;
it = 0;
r(:) = 0;
rt(:) = 0;
rp(:) = 0;
rt2(:) = 0;
rp2(:) = 0;
rtp(:) = 0;
rt3(:) = 0;
rp3(:) = 0;
rt2p(:) = 0;
rtp2(:) = 0;

% Calcaulating surface derivatives
% Loop over harmonics
for n = 0:pxmn
    ih = ih+1;
%   This is assuming that the deformed and reference are at the same tht
%   and phi, which they should be for the first time step.
    Ypcur = Yt{ih};
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
%       Current value of harmonic coefficient        
        f = xmn(it);
%       No need to do the calculation if the current f is zero
        if (f == 0); continue; end
%       Current values of harmonic at degree and order
        Ym = squeeze(Ypcur(im,:,:));
%       Theta derivative
        td1 = ThetDer(Yt,ph,n,m,1);
        
%       Distance from zero
        r = r +  f*Ym;
        
%       Get the derivative in theta and add in
        rt = rt +  f*td1;
%       Get the derivative in phi and add in
        rp = rp +  f*1i*m*Ym;
        
%       Second order derivatives (for curvature
%       Second theta derivative 
        td2 = ThetDer(Yt,ph,n,m,2);

%       Second phi derivative
        rp2 = rp2 +  f*-m^2*Ym;

%       Mixed derivative
        rtp = rtp +  f*1i*m*td1;

%       Second theta derivative
        rt2 = rt2 +  f*td2;

%       Third order derivatives (for curvature derivatives)
%       Third phi derivative
        rp3 = rp3 +  f*-1i*m^3*Ym;

%       Mixed derivatives
        rt2p = rt2p +  f* 1i*m  *td2;
        rtp2 = rtp2 +  f*-1 *m^2*td1;

%       Third theta derivative
        rt3 = rt3 +  f*ThetDer(Yt,ph,n,m,3);
    end
end

% These should be real, but sometimes the imaginary part can creep up a bit
if (max(max(max(abs(imag(rt)))))) > 1e-14
    error('Derivatives in theta are imaginary!')
end

rt   = real(rt);
rt2  = real(rt2);
rtp  = real(rtp);
rt3  = real(rt3);
rt2p = real(rt2p);
rtp2 = real(rtp2);

% Cartesian coordinates (all obtained through the chain rule)
x(1,:,:) = r.*sin(tht).*cos(phi);
x(2,:,:) = r.*sin(tht).*sin(phi);
x(3,:,:) = r.*cos(tht);

% Derivatives in Cartesian (obtained via lots of chain rule)
% Theta
dxt(1,:,:) = rt.*sin(th).*cos(ph) + cos(th).*r.*cos(ph);
dxt(2,:,:) = rt.*sin(th).*sin(ph) + cos(th).*r.*sin(ph);
dxt(3,:,:) = rt.*cos(th) - sin(th).*r;

dxt = real(dxt);

% Phi
dxp(1,:,:) = rp.*sin(th).*cos(ph) - sin(th).*r.*sin(ph);
dxp(2,:,:) = rp.*sin(th).*sin(ph) + sin(th).*r.*cos(ph);
dxp(3,:,:) = rp.*cos(th);

% Normal vector (inward)
for i = 1:nt
    for j = 1:np
        nk(:,i,j) = cross(dxt(:,i,j),dxp(:,i,j));
        nk(:,i,j) = -nk(:,i,j)./norm(nk(:,i,j));
        
    end
end

% Jacobian via fundamental forms
E = squeeze(dot(dxt,dxt));
F = squeeze(dot(dxt,dxp));
G = squeeze(dot(dxp,dxp));
J = sqrt(E.*G - F.^2);

% Second order derivatives
% Second phi
dxp2(1,:,:) = sin(th).*cos(ph).*(rp2 - r) - 2*rp.*sin(th).*sin(ph);
dxp2(2,:,:) = sin(th).*sin(ph).*(rp2 - r) + 2*rp.*sin(th).*cos(ph);
dxp2(3,:,:) = rp2.*cos(th);

% Second theta
dxt2(1,:,:) = sin(th).*cos(ph).*(rt2 - r) + 2*rt.*cos(th).*cos(ph);
dxt2(2,:,:) = sin(th).*sin(ph).*(rt2 - r) + 2*rt.*cos(th).*sin(ph);
dxt2(3,:,:) = cos(th).*(rt2-r) - 2*rt.*sin(th);

% Cross deriv
dxtp(1,:,:) = rtp.*sin(th).*cos(ph) + rp.*cos(th).*cos(ph) - rt.*sin(th).*sin(ph) - r.*cos(th).*sin(ph);
dxtp(2,:,:) = rtp.*sin(th).*sin(ph) + rp.*cos(th).*sin(ph) + rt.*sin(th).*cos(ph) + r.*cos(th).*cos(ph);
dxtp(3,:,:) = rtp.*cos(th) - rp.*sin(th);

% Calculate mean curvature with some fundamental forms (appendices in
% Veeranpani paper)
L = squeeze(dot(dxt2,-nk));
M = squeeze(dot(dxtp,-nk));
N = squeeze(dot(dxp2,-nk));
k = (E.*N - 2*F.*M + G.*L)./(2*J.^2);

% Third Derivatives
% Third ph
dxp3(1,:,:) = sin(th).*cos(ph).*(rp3 - 3*rp) + sin(th).*sin(ph).*(r - 3*rp2);
dxp3(2,:,:) = sin(th).*sin(ph).*(rp3 - 3*rp) + sin(th).*cos(ph).*(3*rp2 - r);
dxp3(3,:,:) = rp3.*cos(th);

% Third theta
dxt3(1,:,:) = sin(th).*cos(ph).*(rt3 - 3*rt) + cos(th).*cos(ph).*(3*rt2 - r);
dxt3(2,:,:) = sin(th).*sin(ph).*(rt3 - 3*rt) + cos(th).*sin(ph).*(3*rt2 - r);
dxt3(3,:,:) = cos(th).*(rt3 - 3*rt) + sin(th).*(r - 3*rt2);

% 2 ph 1 theta
dxtp2(1,:,:)= sin(th).*cos(ph).*(rtp2 - rt) + cos(th).*cos(ph).*(rp2 - r) ... 
            - 2*rtp.*sin(th).*sin(ph) - 2*rp.*cos(th).*sin(ph);
dxtp2(2,:,:)= sin(th).*sin(ph).*(rtp2 - rt) + cos(th).*sin(ph).*(rp2 - r) ... 
            + 2*rtp.*sin(th).*cos(ph) + 2*rp.*cos(th).*cos(ph);
dxtp2(3,:,:)= rtp2.*cos(th) - rp2.*sin(th);

% 1 ph 2 theta
dxt2p(1,:,:)= sin(th).*cos(ph).*(rt2p - rp) + sin(th).*sin(ph).*(r - rt2) ...
            + 2*rtp.*cos(th).*cos(ph) - 2*rt.*cos(th).*sin(ph);
dxt2p(2,:,:)= sin(th).*sin(ph).*(rt2p - rp) + sin(th).*cos(ph).*(rt2 - r) ...
            + 2*rtp.*cos(th).*sin(ph) + 2*rt.*cos(th).*cos(ph);
dxt2p(3,:,:)= cos(th).*(rt2p - rp) - 2*rtp.*sin(th);


for i = 1:nt
    for j = 1:np
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
        dmab(2,1) = (kp - kdR(1,i,j))*gn(2,1) + (k(i,j) - kR(i,j))*dgp(2,1);
        dmab(2,2) = (kp - kdR(1,i,j))*gn(2,2) + (k(i,j) - kR(i,j))*dgp(2,2);
        
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
            + 0.5*es(1)*es(2)*(C*I2-B)*P;
        
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
Yp = SpHarmT(p,th,ph);
f11 = SpT(Yp,tau11,th,ph);
f12 = SpT(Yp,tau12,th,ph);
f21 = SpT(Yp,tau21,th,ph);
f22 = SpT(Yp,tau22,th,ph);
fq1 = SpT(Yp,q1,th,ph);
fq2 = SpT(Yp,q2,th,ph);

for i = 1:nt
    for j =1:np
%       Calculate the partials of tau and q        
        dtau(:) = 0;
        dq(:) = 0;
        ih = 0;
        it = 0;
        for k = 0:p
            ih = ih+1; 
            Ypcur = Yp{ih};
            im = 0;
            for l = -k:k
                it = it+1;
                im = im+1;
%               Current values of harmonic at degree and order
                Ym = squeeze(Ypcur(im,:,:));
                
%               Theta derivative
                td1 = ThetDer(Yp,ph,k,l,1);
                dtau(1,1) = dtau(1,1) + td1(i,j)*f11(it);
                dtau(1,2) = dtau(1,2) + td1(i,j)*f12(it);
                dq(1) = dq(1) + td1(i,j)*fq1(it);
                
%               Phi derivative
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

%% Calculation of velocity constants via fluid problem

% Harmonics at north pole
Ytcr= SpHarmT(pxmn,0,0);

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
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        
%       Velocity at colloc point
        Uc = U + dU*x(:,i,j);        
        
%       Rotation Matrix
        t1 = [cos(phi(j)),-sin(phi(j)),0;sin(phi(j)),cos(phi(j)),0;0,0,1];
        t2 = [cos(-tht(i)),0,sin(-tht(i));0,1,0;-sin(-tht(i)),0,cos(-tht(i))];
        t3 = [cos(-phi(j)),-sin(-phi(j)),0;sin(-phi(j)),cos(-phi(j)),0;0,0,1]; 
        Tx = t1*t2*t3;

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
            end
        end
        
%       Rotated harmonic coefficients of geometry
%       Perhaps we don't need to do this? The Jacobian should be rotation
%       invariant and we rotate the normal vector back anyways?????????????
        xmnr = SpHRot(xmn,phi(j),-tht(i),-phi(j));
        
%       Gauss pts/Normal vector/Jacobian for rotated reference frame
        [xcg, nkg, Jg] = SpHDer(xmnr,Yt,tht,phi);
        
        xcr = [0,0,SpHReconst(xmnr,Ytcr)]';
        
%       Loop over harmonics
        for n = 0:p
%           All spherical harmonics of order n evaluated at integ points
            Ypcur = SpHarm(n,thet,phit);
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
%                         v = Gij(r);
                        v = Tij(r,Tx'*nkg(:,ig,jg));
%                       Integral addition                        
                        At = At + v*Y(ig,jg)*Jg(ig,jg)*ws(ig)*dphi*(1-lam)/(1+lam);
%                       Plain spherical harmonics addition
                        At = At - Y(ig,jg)*4*pi;
                        
%                       Only need to calc B once per colloc point
                        if(n==0)
%                             v = Tij(r,Tx'*nkg(:,ig,jg));
                            v = Gij(r);
%                             bt = bt + v*(U + dU*(Tx'*xcg(:,ig,jg)))*Jg(ig,jg)*ws(ig);
                            bt = bt + v*myf(:,i,j)*Jg(ig,jg)*ws(ig);
                        end
                    end
                end
 
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
                            bt = bt + b(3*ii-2:3*ii)*conj(Y(i,j))*wg(i)*dphi; %!!!!!!!!!!!!!! is this the right thing to be conjugate of
                        end
                    end
                end
                
                A2(row:row+2,col:col+2) = At;
            end
        end
        b2(row:row+2) = bt;
    end
end

[Us,S,V] = svd(A2);
Si = 1./S;
Si(Si>1e10)=0;
Si(end) = 0;
Ai = V*Si*Us';
ut =Ai*b2;

%% Convection of points and new surface constants
% Including getting new angles for convected points, then go back to the
% top

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
ua = -real(ua); %!!!!!!!!!!!!!! Why is this negative?
% xmn = xmn - ua; !!!!!!!!!! This will be money later!
x = x + ua*dt;
r = squeeze(sqrt(x(1,:,:).^2 + x(2,:,:).^2 +x(3,:,:).^2));

xmn = SpT(Yp,r,th,ph)';

clf;
tmpt = linspace(0,pi,100);tmpp = linspace(0,2*pi,100);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmT(p,tmpth,tmpph);
try1 = real(SpHReconst(xmn,Yr));
[Xm,Ym,Zm] = sph2cart(tmpph, tmpth-pi/2, try1);
surf(Xm,Ym,Zm,'edgecolor','none')
axis([-2,2,-2,2,-2,2])
pbaspect([1,1,1])
pause(.01);

end






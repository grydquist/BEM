% !!!!!!! CHECK TO SEE RBC RETURNS TO OG SHAPE
%% Given a deformed and reference configuration, calculate the traction on
% the surface of the object, for an infinitesimal 2d membrane of spherical
% harmonic topology.

% Start with some matrial constants.
B = 0.005;
C = 100;
Eb = 1;%0.0669; % Non-dimensional

%% Defining the shape

% Our grid of points along the surface to get the forces from
    % In the coupled formulation, this grid will be replaced by an initial
    % grid (likely still the Gauss points to start), which will get
    % advected as the cell moves. These points wil then be use to
    % interpolate the forces to the Gauss points for the fluid.
p = 1;
gf = 1;
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);

% Spherical harmonic coefficients defining the shape of the REFERENCE
xmnR = 2*sqrt(pi);
% Spherical harmonic coefficients defining the shape of the DEFORMED
xmn  = 4*sqrt(pi);

% Order of geometric parameterization
pxmn = sqrt(length(xmnR))-1;
% Check that this is an integer
if floor(pxmn) ~= pxmn
    error('Geometric parameterization has the wrong number of constants')
end

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmT(p,th,ph);

%% Getting the surface derivatives, which are needed in a number of places 
% We will need to do this every step for the deformed state, but it only
% needs to be done once for the reference state to get the curvature and
% it's derivatives in theta/phi

% Absolute distance from zero of a given point
r = zeros(nt,np);
rR = zeros(nt,np);

% Derivatives of just harmonics at phi/tht
rt = zeros(nt,np);
rp = rt;

rtR= rt;
rpR= rt;

%   Higher order derivatives
% 2nd
rt2 = rt;
rp2 = rt;
rtp = rt;

rt2R= rt;
rp2R= rt;
rtpR= rt;

% 3rd
rt3 = rt;
rp3 = rt;
rt2p= rt;
rtp2= rt;

rt3R = rt;
rp3R= rt;
rt2pR= rt;
rtp2R= rt;

% Total derivatives, in Cartesian
dxp = zeros(3,nt,np);
dxt = dxp;

dxpR= dxp;
dxtR= dxp;

% Second order
dxp2= dxp;
dxt2= dxp;
dxtp= dxp;

dxp2R= dxp;
dxt2R= dxp;
dxtpR= dxp;

% Third order
dxp3 = dxp;
dxt3 = dxp;
dxt2p= dxp;
dxtp2= dxp;

dxp3R= dxp;
dxt3R= dxp;
dxt2pR= dxp;
dxtp2R= dxp;

ih = 0;
it = 0;

% Now let's actually calculate them
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
        fR= xmnR(it);
%       No need to do the calculation if the current f is zero
        if (f == 0 && fR == 0); continue; end
%       Current values of harmonic at degree and order
        Ym = squeeze(Ypcur(im,:,:));
%       Theta derivative
        td1 = ThetDer(Yt,ph,n,m,1);
        
%       Distance from zero
        r = r +  f*Ym;
        rR= rR+ fR*Ym;
        
%       Get the derivative in theta and add in
        rt = rt +  f*td1;
        rtR= rtR+ fR*td1;
%       Get the derivative in phi and add in
        rp = rp +  f*1i*m*Ym;
        rpR= rpR+ fR*1i*m*Ym;
        
%       Second order derivatives (for curvature
%       Second theta derivative 
        td2 = ThetDer(Yt,ph,n,m,2);

%       Second phi derivative
        rp2 = rp2 +  f*-m^2*Ym;
        rp2R= rp2R+ fR*-m^2*Ym;

%       Mixed derivative
        rtp = rtp +  f*1i*m*td1;
        rtpR= rtpR+ fR*1i*m*td1;

%       Second theta derivative
        rt2 = rt2 +  f*td2;
        rt2R= rt2R+ fR*td2;

%       Third order derivatives (for curvature derivatives)
%       Third phi derivative
        rp3 = rp3 +  f*-1i*m^3*Ym;
        rp3R= rp3R+ fR*-1i*m^3*Ym;

%       Mixed derivatives
        rt2p = rt2p +  f* 1i*m  *td2;
        rtp2 = rtp2 +  f*-1 *m^2*td1;
        rt2pR= rt2pR+ fR* 1i*m  *td2;
        rtp2R= rtp2R+ fR*-1 *m^2*td1;

%       Third theta derivative
        rt3 = rt3 +  f*ThetDer(Yt,ph,n,m,3);
        rt3R= rt3R+ fR*ThetDer(Yt,ph,n,m,3); 
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

rtR   = real(rtR);
rt2R  = real(rt2R);
rtpR  = real(rtpR);
rt3R  = real(rt3R);
rt2pR = real(rt2pR);
rtp2R = real(rtp2R);

% Derivatives in Cartesian (obtained via lots of chain rule)
% Theta
dxt(1,:,:) = rt.*sin(th).*cos(ph) + cos(th).*r.*cos(ph);
dxt(2,:,:) = rt.*sin(th).*sin(ph) + cos(th).*r.*sin(ph);
dxt(3,:,:) = rt.*cos(th) - sin(th).*r;

dxtR(1,:,:) = rtR.*sin(th).*cos(ph) + cos(th).*rR.*cos(ph);
dxtR(2,:,:) = rtR.*sin(th).*sin(ph) + cos(th).*rR.*sin(ph);
dxtR(3,:,:) = rtR.*cos(th) - sin(th).*rR;

dxt = real(dxt);
dxtR = real(dxtR);

% Phi
dxp(1,:,:) = rp.*sin(th).*cos(ph) - sin(th).*r.*sin(ph);
dxp(2,:,:) = rp.*sin(th).*sin(ph) + sin(th).*r.*cos(ph);
dxp(3,:,:) = rp.*cos(th);

dxpR(1,:,:) = rpR.*sin(th).*cos(ph) - sin(th).*rR.*sin(ph);
dxpR(2,:,:) = rpR.*sin(th).*sin(ph) + sin(th).*rR.*cos(ph);
dxpR(3,:,:) = rpR.*cos(th);

% Normal vector (inward)
n = zeros(3,nt,np);
nR = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        n(:,i,j) = cross(dxt(:,i,j),dxp(:,i,j));
        n(:,i,j) = -n(:,i,j)./norm(n(:,i,j));
        
        nR(:,i,j) = cross(dxtR(:,i,j),dxpR(:,i,j));
        nR(:,i,j) = -nR(:,i,j)./norm(nR(:,i,j));
    end
end

% Jacobian via fundamental forms
E = squeeze(dot(dxt,dxt));
F = squeeze(dot(dxt,dxp));
G = squeeze(dot(dxp,dxp));
J = sqrt(E.*G - F.^2);

ER = squeeze(dot(dxtR,dxtR));
FR = squeeze(dot(dxtR,dxpR));
GR = squeeze(dot(dxpR,dxpR));
JR = sqrt(ER.*GR - FR.^2);

% Second order derivatives
% Second phi
dxp2(1,:,:) = sin(th).*cos(ph).*(rp2 - r) - 2*rp.*sin(th).*sin(ph);
dxp2(2,:,:) = sin(th).*sin(ph).*(rp2 - r) + 2*rp.*sin(th).*cos(ph);
dxp2(3,:,:) = rp2.*cos(th);

dxp2R(1,:,:) = sin(th).*cos(ph).*(rp2R - rR) - 2*rpR.*sin(th).*sin(ph);
dxp2R(2,:,:) = sin(th).*sin(ph).*(rp2R - rR) + 2*rpR.*sin(th).*cos(ph);
dxp2R(3,:,:) = rp2R.*cos(th);

% Second theta
dxt2(1,:,:) = sin(th).*cos(ph).*(rt2 - r) + 2*rt.*cos(th).*cos(ph);
dxt2(2,:,:) = sin(th).*sin(ph).*(rt2 - r) + 2*rt.*cos(th).*sin(ph);
dxt2(3,:,:) = cos(th).*(rt2-r) - 2*rt.*sin(th);

dxt2R(1,:,:) = sin(th).*cos(ph).*(rt2R - rR) + 2*rtR.*cos(th).*cos(ph);
dxt2R(2,:,:) = sin(th).*sin(ph).*(rt2R - rR) + 2*rtR.*cos(th).*sin(ph);
dxt2R(3,:,:) = cos(th).*(rt2R-rR) - 2*rtR.*sin(th);

% Cross deriv
dxtp(1,:,:) = rtp.*sin(th).*cos(ph) + rp.*cos(th).*cos(ph) - rt.*sin(th).*sin(ph) - r.*cos(th).*sin(ph);
dxtp(2,:,:) = rtp.*sin(th).*sin(ph) + rp.*cos(th).*sin(ph) + rt.*sin(th).*cos(ph) + r.*cos(th).*cos(ph);
dxtp(3,:,:) = rtp.*cos(th) - rp.*sin(th);

dxtpR(1,:,:) = rtpR.*sin(th).*cos(ph) + rpR.*cos(th).*cos(ph) - rtR.*sin(th).*sin(ph) - rR.*cos(th).*sin(ph);
dxtpR(2,:,:) = rtpR.*sin(th).*sin(ph) + rpR.*cos(th).*sin(ph) + rtR.*sin(th).*cos(ph) + rR.*cos(th).*cos(ph);
dxtpR(3,:,:) = rtpR.*cos(th) - rpR.*sin(th);

% Calculate mean curvature with some fundamental forms (appendices in
% Veeranpani paper)
L = squeeze(dot(dxt2,-n));
M = squeeze(dot(dxtp,-n));
N = squeeze(dot(dxp2,-n));
k = (E.*N - 2*F.*M + G.*L)./(2*J.^2);

LR= squeeze(dot(dxt2R,-nR));
MR= squeeze(dot(dxtpR,-nR));
NR= squeeze(dot(dxp2R,-nR));
kR= (ER.*NR - 2*FR.*MR + GR.*LR)./(2*JR.^2);

% Third Derivatives
% Third ph
dxp3(1,:,:) = sin(th).*cos(ph).*(rp3 - 3*rp) + sin(th).*sin(ph).*(r - 3*rp2);
dxp3(2,:,:) = sin(th).*sin(ph).*(rp3 - 3*rp) + sin(th).*cos(ph).*(3*rp2 - r);
dxp3(3,:,:) = rp3.*cos(th);

dxp3R(1,:,:) = sin(th).*cos(ph).*(rp3R - 3*rpR) + sin(th).*sin(ph).*(rR - 3*rp2R);
dxp3R(2,:,:) = sin(th).*sin(ph).*(rp3R - 3*rpR) + sin(th).*cos(ph).*(3*rp2R - rR);
dxp3R(3,:,:) = rp3R.*cos(th);

% Third theta
dxt3(1,:,:) = sin(th).*cos(ph).*(rt3 - 3*rt) + cos(th).*cos(ph).*(3*rt2 - r);
dxt3(2,:,:) = sin(th).*sin(ph).*(rt3 - 3*rt) + cos(th).*sin(ph).*(3*rt2 - r);
dxt3(3,:,:) = cos(th).*(rt3 - 3*rt) + sin(th).*(r - 3*rt2);

dxt3R(1,:,:) = sin(th).*cos(ph).*(rt3R - 3*rtR) + cos(th).*cos(ph).*(3*rt2R - rR);
dxt3R(2,:,:) = sin(th).*sin(ph).*(rt3R - 3*rtR) + cos(th).*sin(ph).*(3*rt2R - rR);
dxt3R(3,:,:) = cos(th).*(rt3R - 3*rtR) + sin(th).*(rR - 3*rt2R);

% 2 ph 1 theta
dxtp2(1,:,:)= sin(th).*cos(ph).*(rtp2 - rt) + cos(th).*cos(ph).*(rp2 - r) ... 
            - 2*rtp.*sin(th).*sin(ph) - 2*rp.*cos(th).*sin(ph);
dxtp2(2,:,:)= sin(th).*sin(ph).*(rtp2 - rt) + cos(th).*sin(ph).*(rp2 - r) ... 
            + 2*rtp.*sin(th).*cos(ph) + 2*rp.*cos(th).*cos(ph);
dxtp2(3,:,:)= rtp2.*cos(th) - rp2.*sin(th);


dxtp2R(1,:,:)= sin(th).*cos(ph).*(rtp2R - rtR) + cos(th).*cos(ph).*(rp2R - rR) ... 
            - 2*rtpR.*sin(th).*sin(ph) - 2*rpR.*cos(th).*sin(ph);
dxtp2R(2,:,:)= sin(th).*sin(ph).*(rtp2R - rtR) + cos(th).*sin(ph).*(rp2R - rR) ... 
            + 2*rtpR.*sin(th).*cos(ph) + 2*rpR.*cos(th).*cos(ph);
dxtp2R(3,:,:)= rtp2R.*cos(th) - rp2R.*sin(th);

% 1 ph 2 theta
dxt2p(1,:,:)= sin(th).*cos(ph).*(rt2p - rp) + sin(th).*sin(ph).*(r - rt2) ...
            + 2*rtp.*cos(th).*cos(ph) - 2*rt.*cos(th).*sin(ph);
dxt2p(2,:,:)= sin(th).*sin(ph).*(rt2p - rp) + sin(th).*cos(ph).*(rt2 - r) ...
            + 2*rtp.*cos(th).*sin(ph) + 2*rt.*cos(th).*cos(ph);
dxt2p(3,:,:)= cos(th).*(rt2p - rp) - 2*rtp.*sin(th);


dxt2pR(1,:,:)= sin(th).*cos(ph).*(rt2pR - rpR) + sin(th).*sin(ph).*(rR - rt2R) ...
            + 2*rtpR.*cos(th).*cos(ph) - 2*rtR.*cos(th).*sin(ph);
dxt2pR(2,:,:)= sin(th).*sin(ph).*(rt2pR - rpR) + sin(th).*cos(ph).*(rt2R - rR) ...
            + 2*rtpR.*cos(th).*sin(ph) + 2*rtR.*cos(th).*cos(ph);
dxt2pR(3,:,:)= cos(th).*(rt2pR - rpR) - 2*rtpR.*sin(th);
    
    
    
%   Website could be useful, but I couldn't get it working for theta:
%   https://ciks.cbt.nist.gov/~garbocz/paper134/node13.html

%% Calculate the rest of the surface info and then stresses!

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


for i = 1:nt
    for j = 1:np
%       First let's do all the cross products we need
        ct2_p = cross(dxt2(:,i,j),dxp(:,i,j));
        ct_tp = cross(dxt (:,i,j),dxtp(:,i,j));
        ctp_p = cross(dxtp(:,i,j),dxp(:,i,j));
        ct_p2 = cross(dxt (:,i,j),dxp2(:,i,j));
        
        ct2_pR = cross(dxt2R(:,i,j),dxpR(:,i,j));
        ct_tpR = cross(dxtR (:,i,j),dxtpR(:,i,j));
        ctp_pR = cross(dxtpR(:,i,j),dxpR(:,i,j));
        ct_p2R = cross(dxtR (:,i,j),dxp2R(:,i,j));
        
%       Jacobian Partials        
        Jt = dot(-n(:,i,j),ct2_p + ct_tp);
        Jp = dot(-n(:,i,j),ctp_p + ct_p2);
        
        JtR = dot(-nR(:,i,j),ct2_pR + ct_tpR);
        JpR = dot(-nR(:,i,j),ctp_pR + ct_p2R);

%       Normal vector partials
        dnt = (1/J(i,j))*(ct2_p + ct_tp - Jt*-n(:,i,j));
        dnp = (1/J(i,j))*(ctp_p + ct_p2 - Jp*-n(:,i,j));
        
        dntR = (1/JR(i,j))*(ct2_pR + ct_tpR - JtR*-nR(:,i,j));
        dnpR = (1/JR(i,j))*(ctp_pR + ct_p2R - JpR*-nR(:,i,j));

%       Fundamental form partials        
        Et = 2*dot(dxt2(:,i,j),dxt(:,i,j));
        Ep = 2*dot(dxtp(:,i,j),dxt(:,i,j));
        Gt = 2*dot(dxtp(:,i,j),dxp(:,i,j));
        Gp = 2*dot(dxp (:,i,j),dxp2(:,i,j));
        Ft = dot(dxt2(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxtp(:,i,j));
        Fp = dot(dxtp(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxp2(:,i,j));
        Lt = dot(dxt3(:,i,j),-n(:,i,j))  + dot(dxt2(:,i,j),dnt);
        Lp = dot(dxt2p(:,i,j),-n(:,i,j)) + dot(dxt2(:,i,j),dnp);
        Nt = dot(dxtp2(:,i,j),-n(:,i,j)) + dot(dxp2(:,i,j),dnt);
        Np = dot(dxp3(:,i,j),-n(:,i,j))  + dot(dxp2(:,i,j),dnp);
        Mt = dot(dxt2p(:,i,j),-n(:,i,j)) + dot(dxtp(:,i,j),dnt);
        Mp = dot(dxtp2(:,i,j),-n(:,i,j)) + dot(dxtp(:,i,j),dnp);
        
        EtR = 2*dot(dxt2R(:,i,j),dxtR(:,i,j));
        EpR = 2*dot(dxtpR(:,i,j),dxtR(:,i,j));
        GtR = 2*dot(dxtpR(:,i,j),dxpR(:,i,j));
        GpR = 2*dot(dxpR (:,i,j),dxp2R(:,i,j));
        FtR = dot(dxt2R(:,i,j),dxpR(:,i,j)) + dot(dxtR(:,i,j),dxtpR(:,i,j));
        FpR = dot(dxtpR(:,i,j),dxpR(:,i,j)) + dot(dxtR(:,i,j),dxp2R(:,i,j));
        LtR = dot(dxt3R(:,i,j),-nR(:,i,j))  + dot(dxt2R(:,i,j),dntR);
        LpR = dot(dxt2pR(:,i,j),-nR(:,i,j)) + dot(dxt2R(:,i,j),dnpR);
        NtR = dot(dxtp2R(:,i,j),-nR(:,i,j)) + dot(dxp2R(:,i,j),dntR);
        NpR = dot(dxp3R(:,i,j),-nR(:,i,j))  + dot(dxp2R(:,i,j),dnpR);
        MtR = dot(dxt2pR(:,i,j),-nR(:,i,j)) + dot(dxtpR(:,i,j),dntR);
        MpR = dot(dxtp2R(:,i,j),-nR(:,i,j)) + dot(dxtpR(:,i,j),dnpR);

%       Curvature derivatives        
        kt = -2*Jt./J(i,j).*k(i,j) + 0.5*J(i,j).^-2. ... 
                * (Et.*N(i,j) + Nt.*E(i,j) - 2*Ft.*M(i,j) - 2*Mt.*F(i,j) + Gt.*L(i,j) + Lt.*G(i,j));
        kp = -2*Jp./J(i,j).*k(i,j) + 0.5*J(i,j).^-2. ... 
                * (Ep.*N(i,j) + Np.*E(i,j) - 2*Fp.*M(i,j) - 2*Mp.*F(i,j) + Gp.*L(i,j) + Lp.*G(i,j));
            
        ktR = -2*JtR./JR(i,j).*kR(i,j) + 0.5*JR(i,j).^-2. ... 
                * (EtR.*NR(i,j) + NtR.*ER(i,j) - 2*FtR.*MR(i,j) - 2*MtR.*FR(i,j) + GtR.*LR(i,j) + LtR.*GR(i,j));
        kpR = -2*JpR./JR(i,j).*kR(i,j) + 0.5*JR(i,j).^-2. ... 
                * (EpR.*NR(i,j) + NpR.*ER(i,j) - 2*FpR.*MR(i,j) - 2*MpR.*FR(i,j) + GpR.*LR(i,j) + LpR.*GR(i,j));
            

            
%%      Shear transverse
%       Contravariant basis vectors        
        c1 = cross(dxp(:,i,j),-n(:,i,j))/J(i,j);
        c2 = cross(-n(:,i,j),dxt(:,i,j))/J(i,j);
        
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
        dmab(1,1) = (kt - ktR)*gn(1,1) + (k(i,j) - kR(i,j))*dgt(1,1);
        dmab(1,2) = (kt - ktR)*gn(1,2) + (k(i,j) - kR(i,j))*dgt(1,2);
        dmab(2,1) = (kp - kpR)*gn(2,1) + (k(i,j) - kR(i,j))*dgp(2,1);
        dmab(2,2) = (kp - kpR)*gn(2,2) + (k(i,j) - kR(i,j))*dgp(2,2);
        
%       Transverse Shears (Covariant divergence of moment tensor)
        q1(i,j) = Eb*(dmab(1,1) + dmab(2,1) + mab(1,1)*(2*c111(i,j) + c221(i,j)) ...
                + mab(2,1)*(2*c112(i,j) + c222(i,j)) + mab(1,2)*c112(i,j) + mab(2,2)*c122(i,j));
        q2(i,j) = Eb*(dmab(1,2) + dmab(2,2) + mab(1,2)*(c111(i,j) + 2*c221(i,j)) ...
                + mab(2,2)*(c112(i,j) + 2*c222(i,j)) + mab(1,1)*c211(i,j) + mab(2,1)*c221(i,j));
            
%%      In plane Tension
%       Contravariant form of reference tangent
        c1R = cross(dxpR(:,i,j),-nR(:,i,j))/JR(i,j);
        c2R = cross(-nR(:,i,j),dxtR(:,i,j))/JR(i,j);
        
%       Deformation gradient tensor
        Fd = dxt(:,i,j)*c1R' + dxp(:,i,j)*c2R';
        
%       Surface projection operator
        P = eye(3) - n(:,i,j)*n(:,i,j)';
        
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

% Note: we don't really need to save all of the reference states, mostly
% just curvature, curvature derivatives, and reference tangents

%% Do the Spherical Harmonic transform so we can take a surface divergence to get the surface stresses

% Spherical harmonic coefficients for each component of the stress
Yp = SpHarmT(p,th,ph);
f11 = SpT(Yp,tau11,th,ph);
f12 = SpT(Yp,tau12,th,ph);
f21 = SpT(Yp,tau21,th,ph);
f22 = SpT(Yp,tau22,th,ph);
fq1 = SpT(Yp,q1,th,ph);
fq2 = SpT(Yp,q2,th,ph);

dtau  = zeros(2,2);
dq = zeros(2,1);
bm = dtau;
bn = bm;

f = zeros(3,nt,np);
fab = zeros(3,1);
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
        f(:,i,j) = fab(1)*dxt(:,i,j) + fab(2)*dxp(:,i,j) + fab(3)*-n(:,i,j);
    end
end

% Some imaginary components may slip in, so let's make them real
f = real(f);







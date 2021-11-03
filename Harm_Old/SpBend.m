function q = SpBend(fmn, tht, phi, Y, dxt, dxp, kR, kdR)
% Given a reference and deformed state (via SpHarm constants), calculates
% the transverse shear via divergence of the bending moment. Bending moment
% is calculated via the curvature.


% Order of geometric parameterization
pmn = sqrt(length(fmn)) - 1;
Ypcur = Y{1};
[nt, np] = size(squeeze(Ypcur));

% Absolute distance from zero of a given point
r = zeros(nt,np);

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


% Second order
dxp2= zeros(3,nt,np);
dxt2= dxp2;
dxtp= dxp2;

% Third order
dxp3 = dxp2;
dxt3 = dxp2;
dxt2p= dxp2;
dxtp2= dxp2;

ih = 0;
it = 0;

% Loop over harmonics
for n = 0:pmn
    ih = ih+1; 
    Ypcur = Y{ih};
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
%       Current value of harmonic coefficient        
        f = fmn(it);
%       No need to do the calculation if the current f is zero
        if (f == 0); continue; end
%       Current values of harmonic at degree and order
        Ym = squeeze(Ypcur(im,:,:));
%       Theta derivative
        td1 = ThetDer(Y,phi,n,m,1);
        
%       Distance from zero
        r = r + f*Ym;
        
%       Get the derivative in theta and add in
        rt = rt + f*td1;
%       Get the derivative in phi and add in
        rp = rp + f*1i*m*Ym;
        
%       Second order derivatives (for curvature
%       Second theta derivative 
        td2 = ThetDer(Y,phi,n,m,2);

%       Second phi derivative
        rp2 = rp2 + f*-m^2*Ym;

%       Mixed derivative
        rtp = rtp + f*1i*m*td1;

%       Second theta derivative
        rt2 = rt2 + f*td2;
        
%       Third order derivatives (for curvature derivatives)
%       Third phi derivative
        rp3 = rp3 + f*-1i*m^3*Ym;

%       Mixed derivatives
        rt2p = rt2p + f* 1i*m  *td2;
        rtp2 = rtp2 + f*-1 *m^2*td1;

%       Third theta derivative
        rt3 = rt3 + f*ThetDer(Y,phi,n,m,3);
        
    end
end

if (max(max(max(abs(imag(rt)))))) > 1e-14
    error('Derivatives in theta are imaginary!')
end
rt = real(rt);

% Normal vector (inward)
n = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        n(:,i,j) = cross(dxt(:,i,j),dxp(:,i,j));
        n(:,i,j) = -n(:,i,j)./norm(n(:,i,j));
    end
end

% Jacobian via fundamental forms
E = squeeze(dot(dxt,dxt));
F = squeeze(dot(dxt,dxp));
G = squeeze(dot(dxp,dxp));
J = sqrt(E.*G - F.^2);

if (max(max(max(abs(imag(rt2)))))) > 1e-14
    error('Second derivatives in theta are imaginary!')
end

rt2 = real(rt2);
rtp = real(rtp);

% Second phi
dxp2(1,:,:) = sin(tht).*cos(phi).*(rp2 - r) - 2*rp.*sin(tht).*sin(phi);
dxp2(2,:,:) = sin(tht).*sin(phi).*(rp2 - r) + 2*rp.*sin(tht).*cos(phi);
dxp2(3,:,:) = rp2.*cos(tht);

% Second theta
dxt2(1,:,:) = sin(tht).*cos(phi).*(rt2 - r) + 2*rt.*cos(tht).*cos(phi);
dxt2(2,:,:) = sin(tht).*sin(phi).*(rt2 - r) + 2*rt.*cos(tht).*sin(phi);
dxt2(3,:,:) = cos(tht).*(rt2-r) - 2*rt.*sin(tht);

% Cross deriv
dxtp(1,:,:) = rtp.*sin(tht).*cos(phi) + rp.*cos(tht).*cos(phi) - rt.*sin(tht).*sin(phi) - r.*cos(tht).*sin(phi);
dxtp(2,:,:) = rtp.*sin(tht).*sin(phi) + rp.*cos(tht).*sin(phi) + rt.*sin(tht).*cos(phi) + r.*cos(tht).*cos(phi);
dxtp(3,:,:) = rtp.*cos(tht) - rp.*sin(tht);

% Calculate curvature with some fundamental forms (appendices in
% Veeranpani paper)
L = squeeze(dot(dxt2,-n));
M = squeeze(dot(dxtp,-n));
N = squeeze(dot(dxp2,-n));
k = (E.*N - 2*F.*M + G.*L)./(2*J.^2);    

rt3  = real(rt3);
rt2p = real(rt2p);
rtp2 = real(rtp2);

%   Third phi
dxp3(1,:,:) = sin(tht).*cos(phi).*(rp3 - 3*rp) + sin(tht).*sin(phi).*(r - 3*rp2);
dxp3(2,:,:) = sin(tht).*sin(phi).*(rp3 - 3*rp) + sin(tht).*cos(phi).*(3*rp2 - r);
dxp3(3,:,:) = rp3.*cos(tht);

%   Third theta
dxt3(1,:,:) = sin(tht).*cos(phi).*(rt3 - 3*rt) + cos(tht).*cos(phi).*(3*rt2 - r);
dxt3(2,:,:) = sin(tht).*sin(phi).*(rt3 - 3*rt) + cos(tht).*sin(phi).*(3*rt2 - r);
dxt3(3,:,:) = cos(tht).*(rt3 - 3*rt) + sin(tht).*(r - 3*rt2);

%   2 phi 1 theta
dxtp2(1,:,:)= sin(tht).*cos(phi).*(rtp2 - rt) + cos(tht).*cos(phi).*(rp2 - r) ... 
            - 2*rtp.*sin(tht).*sin(phi) - 2*rp.*cos(tht).*sin(phi);
dxtp2(2,:,:)= sin(tht).*sin(phi).*(rtp2 - rt) + cos(tht).*sin(phi).*(rp2 - r) ... 
            + 2*rtp.*sin(tht).*cos(phi) + 2*rp.*cos(tht).*cos(phi);
dxtp2(3,:,:)= rtp2.*cos(tht) - rp2.*sin(tht);

%   1 phi 2 theta
dxt2p(1,:,:)= sin(tht).*cos(phi).*(rt2p - rp) + sin(tht).*sin(phi).*(r - rt2) ...
            + 2*rtp.*cos(tht).*cos(phi) - 2*rt.*cos(tht).*sin(phi);
dxt2p(2,:,:)= sin(tht).*sin(phi).*(rt2p - rp) + sin(tht).*cos(phi).*(rt2 - r) ...
            + 2*rtp.*cos(tht).*sin(phi) + 2*rt.*cos(tht).*cos(phi);
dxt2p(3,:,:)= cos(tht).*(rt2p - rp) - 2*rtp.*sin(tht);

%   Derivatives of normal vectors
dnp = zeros(3,nt,np);
dnt = dnp;

%   Fundamental form derivatives (could be a bit overkill)
Et = zeros(nt,np);
Ep = Et;
Ft = Et;
Fp = Et;
Gt = Et;
Gp = Et;
Lt = Et;
Lp = Et;
Mt = Et;
Mp = Et;
Nt = Et;
Np = Et;
Jt = Et;
Jp = Et;
g = zeros(2,2,nt,np);
dmab = zeros(2,2);
q = zeros(2,nt,np);

for i = 1:nt
    for j = 1:np
%       First let's do all the cross products we need
        ct2_p = cross(dxt2(:,i,j),dxp(:,i,j));
        ct_tp = cross(dxt (:,i,j),dxtp(:,i,j));
        ctp_p = cross(dxtp(:,i,j),dxp(:,i,j));
        ct_p2 = cross(dxt (:,i,j),dxp2(:,i,j));
        
%       Contravariant vectors        
        c1 = cross(dxp(:,i,j),-n(:,i,j))/J(i,j);
        c2 = cross(-n(:,i,j),dxt(:,i,j))/J(i,j);
        
        Jt(i,j) = dot(-n(:,i,j),ct2_p + ct_tp);
        Jp(i,j) = dot(-n(:,i,j),ctp_p + ct_p2);

        dnt(:,i,j) = (1/J(i,j))*(ct2_p + ct_tp ...
                   - Jt(i,j)*-n(:,i,j));
        dnp(:,i,j) = (1/J(i,j))*(ctp_p + ct_p2 ...
                   - Jp(i,j)*-n(:,i,j));

        Et(i,j) = 2*dot(dxt2(:,i,j),dxt(:,i,j));
        Ep(i,j) = 2*dot(dxtp(:,i,j),dxt(:,i,j));
        Gt(i,j) = 2*dot(dxtp(:,i,j),dxp(:,i,j));
        Gp(i,j) = 2*dot(dxp (:,i,j),dxp2(:,i,j));
        Ft(i,j) = dot(dxt2(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxtp(:,i,j));
        Fp(i,j) = dot(dxtp(:,i,j),dxp(:,i,j)) + dot(dxt(:,i,j),dxp2(:,i,j));
        Lt(i,j) = dot(dxt3(:,i,j),-n(:,i,j))  + dot(dxt2(:,i,j),dnt(:,i,j));
        Lp(i,j) = dot(dxt2p(:,i,j),-n(:,i,j)) + dot(dxt2(:,i,j),dnp(:,i,j));
        Nt(i,j) = dot(dxtp2(:,i,j),-n(:,i,j)) + dot(dxp2(:,i,j),dnt(:,i,j));
        Np(i,j) = dot(dxp3(:,i,j),-n(:,i,j))  + dot(dxp2(:,i,j),dnp(:,i,j));
        Mt(i,j) = dot(dxt2p(:,i,j),-n(:,i,j)) + dot(dxtp(:,i,j),dnt(:,i,j));
        Mp(i,j) = dot(dxtp2(:,i,j),-n(:,i,j)) + dot(dxtp(:,i,j),dnp(:,i,j));
        g(:,:,i,j) = [E(i,j),F(i,j);F(i,j),G(i,j)];
%       Contravariant g        
        gn = inv(g(:,:,i,j));
%       Derivatives of CONTRAVARIANT g!        
        dgt = -gn*[Et(i,j),Ft(i,j);Ft(i,j),Gt(i,j)]*gn;
        dgp = -gn*[Ep(i,j),Fp(i,j);Fp(i,j),Gp(i,j)]*gn;
        
%       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
        c111 = dot(dxt2(:,i,j),c1);
        c112 = dot(dxtp(:,i,j),c1);
        c122 = dot(dxp2(:,i,j),c1);
        c222 = dot(dxp2(:,i,j),c2);
        c221 = dot(dxtp(:,i,j),c2);
        c211 = dot(dxt2(:,i,j),c2);
        
%       Curvature derivatives
        kt = -2*Jt(i,j)/J(i,j)*k(i,j) + 0.5*J(i,j)^-2 ...
           * (Et(i,j)*N(i,j) + Nt(i,j)*E(i,j) - 2*Ft(i,j)*M(i,j) ...
           - 2*Mt(i,j)*F(i,j) + Gt(i,j)*L(i,j) + Lt(i,j)*G(i,j));
        kp = -2*Jp(i,j)/J(i,j)*k(i,j) + 0.5*J(i,j)^-2 ...
           * (Ep(i,j)*N(i,j) + Np(i,j)*E(i,j) - 2*Fp(i,j)*M(i,j) ...
           - 2*Mp(i,j)*F(i,j) + Gp(i,j)*L(i,j) + Lp(i,j)*G(i,j));

%       Moments and their needed partials
        mab = (k(i,j) - kR(i,j))*gn;
        dmab(1,1) = (kt - kdR(1,i,j))*gn(1,1) + (k(i,j) - kR(i,j))*dgt(1,1);
        dmab(1,2) = (kt - kdR(1,i,j))*gn(1,2) + (k(i,j) - kR(i,j))*dgt(1,2);
        dmab(2,1) = (kp - kdR(2,i,j))*gn(2,1) + (k(i,j) - kR(i,j))*dgp(2,1);
        dmab(2,2) = (kp - kdR(2,i,j))*gn(2,2) + (k(i,j) - kR(i,j))*dgp(2,2);
        
%       Transverse Shears (Covariant divergence of moment tensor)
%       May not be right!!!!!!!!!!
        q(1,i,j) = dmab(1,1) + dmab(2,1) + mab(1,1)*(2*c111 + c221) ...
                 + mab(2,1)*(2*c112 + c222) + mab(1,2)*c112 + mab(2,2)*c122;
        q(2,i,j) = dmab(1,2) + dmab(2,2) + mab(1,2)*(c111 + 2*c221) ...
                 + mab(2,2)*(c112 + 2*c222) + mab(1,1)*c211 + mab(2,1)*c221;
        
             
             % !!!!!!!!!!!!!!!!!!!!!!! LOTS OF unnecesarily large matrices,
             % e.g. g, dgt, maybe even E, Et, etc.
             
    end
end



end

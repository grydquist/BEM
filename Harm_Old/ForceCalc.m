function f = ForceCalc(fmn, ftot, Y, tht, phi, dxt, dxp)
% Calculates the force given geometry and SpHarm constants of stress tensor

% Info on points
Ypcur = Y{1};
[nt, np] = size(squeeze(Ypcur));

% Unpack constants for forces
f11 = ftot(:,1);
f12 = ftot(:,2);
f21 = ftot(:,3);
f22 = ftot(:,4);
fq1 = ftot(:,5);
fq2 = ftot(:,6);
pxmn = sqrt(length(fmn)) - 1;
p    = sqrt(length(f11))-1;

% Actual values of the forces at the points
tau11 = SpHReconst(f11,Y);
tau12 = SpHReconst(f12,Y);
tau21 = SpHReconst(f21,Y);
tau22 = SpHReconst(f22,Y);
q1 = SpHReconst(fq1,Y);
q2 = SpHReconst(fq2,Y);

% Absolute distance
r = zeros(nt,np);

% Derivatives of just harmonics at phi/tht
rt = zeros(nt,np);
rp = rt;
%   Higher order derivatives
% 2nd
rt2 = rt;
rp2 = rt;
rtp = rt;

% Second order
dxp2= zeros(3,nt,np);
dxt2= dxp2;
dxtp= dxp2;

% Calculate the rest of the derivatives
% Loop over harmonics
ih = 0;
it = 0;
for n = 0:pxmn
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
    end
end

% Normal (inward, for consistency)
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

dtau = zeros(2,2);
bm = dtau;
bn = bm;
dq = zeros(2,1);

f = zeros(3,nt,np);

%ALL OF THIS NEEDS TO BE CHECKED

for i = 1:nt
    for j = 1:np
        
%       Contravariant vectors        
        c1 = cross(dxp(:,i,j),-n(:,i,j))/J(i,j);
        c2 = cross(-n(:,i,j),dxt(:,i,j))/J(i,j);

%       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
        c111 = dot(dxt2(:,i,j),c1);
        c112 = dot(dxtp(:,i,j),c1);
        c122 = dot(dxp2(:,i,j),c1);
        c222 = dot(dxp2(:,i,j),c2);
        c221 = dot(dxtp(:,i,j),c2);
        c211 = dot(dxt2(:,i,j),c2);
        
%       Calculate the partials of tau and q        
        dtau(:) = 0;
        dq(:) = 0;
        ih = 0;
        it = 0;
        for k = 0:p
            ih = ih+1; 
            Ypcur = Y{ih};
            im = 0;
            for l = -k:k
                it = it+1;
                im = im+1;
%               Current values of harmonic at degree and order
                Ym = squeeze(Ypcur(im,:,:));
                
%               Theta derivative !!!!!!!!!! Horribly, horribly inefficient
%               AND ALSO, LOOK AT THETDER AGAIN. This won't work if every
%               grid point is at a unique theta
                td1 = ThetDer(Y,phi,k,l,1);
                dtau(1,1) = dtau(1,1) + td1(i)*f11(it);
                dtau(1,2) = dtau(1,2) + td1(i)*f12(it);
                dq(1) = dq(1) + td1(i)*fq1(it);
                
%               Phi derivative
                dtau(2,1) = dtau(2,1) + f21(it)*1i*l*Ym(i,j);
                dtau(2,2) = dtau(2,2) + f22(it)*1i*l*Ym(i,j);
                dq(2) = dq(2) + fq2(it)*1i*l*Ym(i,j);
            end
        end
%       Covariant metric tensor
        gv = [E(i,j),F(i,j);F(i,j),G(i,j)];

%       Contravariant metric 
        gn = inv(gv);

%       Covariant curvature tensor
        bv = [L(i,j),M(i,j);M(i,j),N(i,j)];
        
%       Mixed curvature tensor (raise the column)
        bm(1,1) = bv(1,1)*gn(1,1) + bv(1,2)*gn(2,1);
        bm(1,2) = bv(1,1)*gn(1,2) + bv(1,2)*gn(2,2);
        bm(2,1) = bv(2,1)*gn(1,1) + bv(2,2)*gn(2,1);
        bm(2,2) = bv(2,1)*gn(1,2) + bv(2,2)*gn(2,2);
        
%       Contravariant curvature (raise the row) !!!!!!!!!! NOT needed, but
%       check anyway
        bn(1,1) = bm(1,1)*gn(1,1) + bm(2,1)*gn(2,1);
        bn(1,2) = bm(1,2)*gn(1,1) + bm(2,1)*gn(2,1);
        bn(2,1) = bm(1,1)*gn(1,2) + bm(2,1)*gn(2,2);
        bn(2,2) = bm(1,2)*gn(1,2) + bm(2,2)*gn(2,2);
        
%       Covariant divergence of tau in theta, then phi
        cvt = dtau(1,1) + dtau(2,1) + tau11(i,j)*(2*c111 + c221) ...
            + tau21(i,j)*(2*c112 + c222) + tau12(i,j)*c112 + tau22(i,j)*c122;
        cvp = dtau(1,2) + dtau(2,2) + tau12(i,j)*(c111 + 2*c221) ...
            + tau22(i,j)*(c112 + 2*c222) + tau11(i,j)*c211 + tau21(i,j)*c221;
%       Covariant divergence of Q        
        cq = dq(1) + dq(2) + c111*q1(i,j) + c112*q2(i,j) + c221*q1(i,j) + c222*q2(i,j);
        
        f(1,i,j) = bm(1,1)*q1(i,j) + bm(2,1)*q2(i,j) - cvt;
        f(2,i,j) = bm(1,2)*q1(i,j) + bm(2,2)*q2(i,j) - cvp;
        f(3,i,j) = -cq - tau11(i,j)*bv(1,1) - tau12(i,j)*bv(1,2) - tau21(i,j)*bv(2,1) - tau22(i,j)*bv(2,2);
    end
end









end
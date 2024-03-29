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
% B = 0.005;
% C = 100;
B=.001;% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C=100;

% Bending modulus
Eb = 0;%1;%0.0669; % Non-dimensional !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Total time steps
NT = 1400;
% dt = 0.00001;
dt = 0.005;

% Velocity and gradient
U = [0;0;0];

shft = [cos(pi/4),0,sin(pi/4);0,1,0;-sin(pi/4),0,cos(pi/4)];
dU = [0,0,1;0,0,0;.0,0,0]/4/pi;
% dU = shft'*[0,0,1;0,0,0;1,0,0]*shft;
%% Reference Shape Definition

% Order of the force and velocity transforms
p = 4;

% I call factorial a lot and it makes the most sense to do all the ones
% I'll need right off the bat and store them (goes up to 2p)
myfacs = zeros(2*p+1,1);
for i = 0:2*p
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

% Spherical harmonic coefficients defining the shape of ref & def:
xmnR = zeros((p+1)^2,1);
xmnR(1) = 2*sqrt(pi);

% Get harmonics for the individual parts of the position vector (would be
% good to do this for the reference too, but don't need to at this point)
rrr = SpHReconst(xmnR,Yt);
x1 = rrr.*sin(th).*cos(ph);
x2 = rrr.*sin(th).*sin(ph);
x3 = rrr.*cos(th);

% x1mn = SpT(Yt,x1,th,ph);
% x2mn = SpT(Yt,x2,th,ph);
% x3mn = SpT(Yt,x3,th,ph);

% For the RBC, fully defined at p = 5
[x1mn,x2mn,x3mn] = RBCcoeffs(Yt,th,ph);
ix3mn = x3mn;

%% Preallocation of all the solid things I need

% Actual locations/areas
x = zeros(3,nt,np);
J = zeros(nt,np);

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

% Moment
dmab = zeros(2,2);

dtau  = zeros(2,2);
dq = zeros(2,1);
bm = dtau;
bn = bm;

myf = zeros(3,nt,np);
Nmyf = myf;
Nmyfab = myf;

% Reference state variables that need saving
c1R = zeros(3,nt,np,nt,np);
c2R = c1R;
dc1tR= c1R;
dc1pR= c1R;
dc2tR= c1R;
dc2pR= c1R;
kR = zeros(nt,np,nt,np);
kdR = zeros(2,nt,np,nt,np);

% Get all the derivatives of the Spherical Harmonics up front
Ytd1 = zeros(nt,np,(p+1)^2);
Ytd2 = Ytd1;
Ytd3 = Ytd1;

it = 0;
im = 0;
for n = 0:p
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
        
%       Theta derivative
        Ytd1(:,:,it) = ThetDer(Yt,ph,n,m,1);
        
%       Second theta derivative 
        Ytd2(:,:,it) = ThetDer(Yt,ph,n,m,2);

%       Third theta derivative
        Ytd3(:,:,it) = ThetDer(Yt,ph,n,m,3);
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
trcvel = zeros(1,NT);
trcvelz = zeros(1,NT);

% Spherical harmonic evaluated at right hand side of sphere
Ytrc = SpHarmT(p,pi/2,0,myfacs);

%% Time stepping loop
for cts = 1:NT
%% Calculation of residual force at convected points
% After calculation, do SpHarm transform
if(cts==25)
%     dU(:)=0;
end %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Get Geometric information for current configuration (Perhaps the best way
% to do this is with a subroutine that takes in the SpHarm consts for the
% def, and infor for the ref, and material constants and just spits out the
% forces. In FORTRAN this information would be packaged up in objects, and
% we could just pass pointers to the objects!)

x(1,:,:) = real(SpHReconst(x1mn,Yt));
x(2,:,:) = real(SpHReconst(x2mn,Yt));
x(3,:,:) = real(SpHReconst(x3mn,Yt));

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

%               Using three Spherical harmonics for the three vector coordinates
%               means that the reconstruction will give e.g. the x value in the
%               NON-ROTATED referece frame, USING THE ROTATED reference frame
%               angles (for example, say we rotate the reference 90 deg so the pole
%               is at the x axis in the non rotated frame. Then, an angle of 0,0 in
%               the rotated frame for the x coord returns the value of X at 90,0).

            end
        end
        
%       The location of the new pole, should just be location of
%       point at current theta/phi

%       Harms on rotated grid
        Ypt = SpHarmT(p,thet,phit,myfacs);

%       Forces on rotated grid.

        ih = 0;
        it = 0;
        dxt(:) = 0;
        dxp(:) = 0;
        dxt2(:) = 0;
        dxp2(:) = 0;
        dxtp(:) = 0;
        dxt3(:) = 0;
        dxp3(:) = 0;
        dxt2p(:) = 0;
        dxtp2(:) = 0;

%       Surface derivatives: loop over harmonics        
        for n = 0:p
            ih = ih+1;
        %   This is assuming that the deformed and reference are at the same tht
        %   and phi, which they should be, as they aren't angles necessarily
            Ypcur = Yt{ih};
            im = 0;
            for m = -n:n
                im = im + 1;
                it = it + 1;
        %       Current value of harmonic coefficient        
                f1 = x1mnr(it);
                f2 = x2mnr(it);
                f3 = x3mnr(it);
        %       Current values of harmonic at degree and order
                Ym = squeeze(Ypcur(im,:,:));
        %       Theta derivative
                td1 = Ytd1(:,:,it);

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
            end
        end
%       These shouldn't be imaginary, but some creeps in        
        dxt   = real(dxt);
        dxt2  = real(dxt2);
        dxt3  = real(dxt3);
        
        dxp   = real(dxp);
        dxp2  = real(dxp2);
        dxp3  = real(dxp3);
        
        dxtp  = real(dxtp);
        dxt2p = real(dxt2p);
        dxtp2 = real(dxtp2);
        
%       Here come the stresses!        
        for i2 = 1:nt
            for j2 = 1:np

        %       Normal vector (inward)
                nk(:,i2,j2) = cross(dxt(:,i2,j2),dxp(:,i2,j2));
                nk(:,i2,j2) = -nk(:,i2,j2)./norm(nk(:,i2,j2));
                nk(:,i2,j2) = real(nk(:,i2,j2));

        %       Jacobian via fundamental forms
                E = dot(dxt(:,i2,j2),dxt(:,i2,j2));
                F = dot(dxt(:,i2,j2),dxp(:,i2,j2));
                G = dot(dxp(:,i2,j2),dxp(:,i2,j2));
                J(i2,j2) = sqrt(E.*G - F.^2);

        %       Calculate mean curvature with some fundamental forms (appendices in
        %       Veeranpani paper)
                L = dot(dxt2(:,i2,j2),-nk(:,i2,j2));
                M = dot(dxtp(:,i2,j2),-nk(:,i2,j2));
                N = dot(dxp2(:,i2,j2),-nk(:,i2,j2));
                k = (E*N - 2*F*M + G*L)/(2*J(i2,j2)^2);

        %       First let's do all the cross products we need
                ct2_p = cross(dxt2(:,i2,j2),dxp(:,i2,j2));
                ct_tp = cross(dxt (:,i2,j2),dxtp(:,i2,j2));
                ctp_p = cross(dxtp(:,i2,j2),dxp(:,i2,j2));
                ct_p2 = cross(dxt (:,i2,j2),dxp2(:,i2,j2));

        %       Jacobian Partials        
                Jt = dot(-nk(:,i2,j2),ct2_p + ct_tp);
                Jp = dot(-nk(:,i2,j2),ctp_p + ct_p2);

        %       Normal vector partials outward(?)
                dnt = (1/J(i2,j2))*(ct2_p + ct_tp - Jt*-nk(:,i2,j2));
                dnp = (1/J(i2,j2))*(ctp_p + ct_p2 - Jp*-nk(:,i2,j2));

        %       Fundamental form partials        
                Et = 2*dot(dxt2(:,i2,j2),dxt(:,i2,j2));
                Ep = 2*dot(dxtp(:,i2,j2),dxt(:,i2,j2));
                Gt = 2*dot(dxtp(:,i2,j2),dxp(:,i2,j2));
                Gp = 2*dot(dxp (:,i2,j2),dxp2(:,i2,j2));
                Ft = dot(dxt2(:,i2,j2),dxp(:,i2,j2)) + dot(dxt(:,i2,j2),dxtp(:,i2,j2));
                Fp = dot(dxtp(:,i2,j2),dxp(:,i2,j2)) + dot(dxt(:,i2,j2),dxp2(:,i2,j2));
                Lt = dot(dxt3(:,i2,j2),-nk(:,i2,j2))  + dot(dxt2(:,i2,j2),dnt);
                Lp = dot(dxt2p(:,i2,j2),-nk(:,i2,j2)) + dot(dxt2(:,i2,j2),dnp);
                Nt = dot(dxtp2(:,i2,j2),-nk(:,i2,j2)) + dot(dxp2(:,i2,j2),dnt);
                Np = dot(dxp3(:,i2,j2),-nk(:,i2,j2))  + dot(dxp2(:,i2,j2),dnp);
                Mt = dot(dxt2p(:,i2,j2),-nk(:,i2,j2)) + dot(dxtp(:,i2,j2),dnt);
                Mp = dot(dxtp2(:,i2,j2),-nk(:,i2,j2)) + dot(dxtp(:,i2,j2),dnp);

        %       Curvature derivatives        
                kt = -2*Jt./J(i2,j2).*k + 0.5*J(i2,j2).^-2. ... 
                        * (Et.*N + Nt.*E - 2*Ft.*M - 2*Mt.*F + Gt.*L + Lt.*G);
                kp = -2*Jp./J(i2,j2).*k + 0.5*J(i2,j2).^-2. ... 
                        * (Ep.*N + Np.*E - 2*Fp.*M - 2*Mp.*F + Gp.*L + Lp.*G);   

        %       Metric tensor
                g = [E,F;F,G];
        %       Contravariant g
                gn = inv(g);
        %       Partials of COTNRAVARIANT metric tensor
                dgt = -gn*[Et,Ft;Ft,Gt]*gn;
                dgp = -gn*[Ep,Fp;Fp,Gp]*gn;    

        %       Shear transverse
        %       Contravariant basis vectors, raise index 
                c1 = dxt(:,i2,j2)*gn(1,1) + dxp(:,i2,j2)*gn(1,2);
                c2 = dxt(:,i2,j2)*gn(2,1) + dxp(:,i2,j2)*gn(2,2);

        %       !!!!!!!!!! HACKY FIRST TIME STEP
                if(cts == 1)
                    c1R(:,i,j,i2,j2) = c1;
                    c2R(:,i,j,i2,j2) = c2;
                    kR(i,j,i2,j2) = k;
                    kdR(1,i,j,i2,j2) = kt;
                    kdR(2,i,j,i2,j2) = kp;
        %           For derivatives of tau, we need dcR/dtht and dcR/dphi
        %           Chain rule across raise index operation
                    dc1tR(:,i,j,i2,j2) = dxt2(:,i2,j2)*gn(1,1) + dxtp(:,i2,j2)*gn(2,1)...
                                 + dxt(:,i2,j2)*dgt(1,1) + dxp(:,i2,j2)*dgt(2,1);

                    dc1pR(:,i,j,i2,j2) = dxtp(:,i2,j2)*gn(1,1) + dxp2(:,i2,j2)*gn(2,1)...
                                 + dxt(:,i2,j2)*dgp(1,1) + dxp(:,i2,j2)*dgp(2,1);

                    dc2tR(:,i,j,i2,j2) = dxt2(:,i2,j2)*gn(1,2) + dxtp(:,i2,j2)*gn(2,2)...
                                 + dxt(:,i2,j2)*dgt(1,2) + dxp(:,i2,j2)*dgt(2,2);

                    dc2pR(:,i,j,i2,j2) = dxtp(:,i2,j2)*gn(1,2) + dxp2(:,i2,j2)*gn(2,2)...
                                 + dxt(:,i2,j2)*dgp(1,2) + dxp(:,i2,j2)*dgp(2,2);
                end

        %       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
                c111 = dot(dxt2(:,i2,j2),c1);
                c112 = dot(dxtp(:,i2,j2),c1);
                c122 = dot(dxp2(:,i2,j2),c1);
                c222 = dot(dxp2(:,i2,j2),c2);
                c221 = dot(dxtp(:,i2,j2),c2);
                c211 = dot(dxt2(:,i2,j2),c2);

        %       Moments and their needed partials
                mab = (k - kR(i,j,i2,j2))*gn;
                dmab(1,1) = (kt - kdR(1,i,j,i2,j2))*gn(1,1) + (k - kR(i,j,i2,j2))*dgt(1,1);
                dmab(1,2) = (kt - kdR(1,i,j,i2,j2))*gn(1,2) + (k - kR(i,j,i2,j2))*dgt(1,2);
                dmab(2,1) = (kp - kdR(2,i,j,i2,j2))*gn(2,1) + (k - kR(i,j,i2,j2))*dgp(2,1);
                dmab(2,2) = (kp - kdR(2,i,j,i2,j2))*gn(2,2) + (k - kR(i,j,i2,j2))*dgp(2,2);

        %       Transverse Shears (Covariant divergence of moment tensor)
                q1 = Eb*(dmab(1,1) + dmab(2,1) + mab(1,1)*(2*c111 + c221) ...
                        + mab(2,1)*(2*c112 + c222) + mab(1,2)*c112 + mab(2,2)*c122);
                q2 = Eb*(dmab(1,2) + dmab(2,2) + mab(1,2)*(c111 + 2*c221) ...
                        + mab(2,2)*(c112 + 2*c222) + mab(1,1)*c211 + mab(2,1)*c221);

        %%      In plane Tension
        %       Deformation gradient tensor
                Fd = dxt(:,i2,j2)*c1R(:,i,j,i2,j2)' + dxp(:,i2,j2)*c2R(:,i,j,i2,j2)';

        %       Surface projection operator
                P = eye(3) - nk(:,i2,j2)*nk(:,i2,j2)';

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
                Fdv = [dxt(:,i2,j2)'*Fd*dxt(:,i2,j2),dxt(:,i2,j2)'*Fd*dxp(:,i2,j2)
                       dxp(:,i2,j2)'*Fd*dxt(:,i2,j2),dxp(:,i2,j2)'*Fd*dxp(:,i2,j2)];
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
                dFdt = dxt2(:,i2,j2)*c1R(:,i,j,i2,j2)'  + dxtp(:,i2,j2)*c2R(:,i,j,i2,j2)' ...
                    + dxt(:,i2,j2)*dc1tR(:,i,j,i2,j2)' + dxp(:,i2,j2)*dc2tR(:,i,j,i2,j2)';

                dFdp = dxtp(:,i2,j2)*c1R(:,i,j,i2,j2)'  + dxp2(:,i2,j2)*c2R(:,i,j,i2,j2)' ...
                    + dxt(:,i2,j2)*dc1pR(:,i,j,i2,j2)' + dxp(:,i2,j2)*dc2pR(:,i,j,i2,j2)';

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
                dPt = -dnt*-nk(:,i2,j2)' + nk(:,i2,j2)*dnt';
                dPp = -dnp*-nk(:,i2,j2)' + nk(:,i2,j2)*dnp';

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

%%      And finally the forces!
                
        %       Metric tensor
                g = [E,F;F,G];
        %       Contravariant g
                gn = inv(g);

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
                fab(1) = bm(1,1)*q1 + bm(2,1)*q2 - cvt;
                fab(2) = bm(1,2)*q1 + bm(2,2)*q2 - cvp;
                fab(3) = -cq - tau11*bv(1,1) - tau12*bv(1,2) - tau21*bv(2,1) - tau22*bv(2,2);

        %       These are in terms of (contravariant) surface vectors, put them into Cartesian
                myf(:,i2,j2) = fab(1)*dxt(:,i2,j2) + fab(2)*dxp(:,i2,j2) + fab(3)*-nk(:,i2,j2);
            end
        end

        % Some imaginary components may slip in, so let's make them real
        Nmyf = real(myf);

        xcr = x(:,i,j);

%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        
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
                        r = (xcr-xcg(:,ig,jg));
                        v = Tij(r,nk(:,ig,jg));
%                       Integral addition                        
                        At = At + v*Y(ig,jg)*J(ig,jg)*ws(ig)*dphi;
                        
%                       Only need to calc B once per colloc point
                        if(n==0)
                            v = Gij(r);
                            bt = bt + v*Nmyf(:,ig,jg)*J(ig,jg)*ws(ig);
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
% A2 = -4*pi*eye(3*(p+1)^2);
ut = A2\b2;

%% Convection of points and new surface constants

% Throw out highest order coefficients (give lots of error, unstable)
ut(3*p^2+1:end) = 0;

u1 = ut(1:3:end);
u2 = ut(2:3:end);
u3 = ut(3:3:end);
x1mn = x1mn + u1*dt;
x2mn = x2mn + u2*dt;
x3mn = x3mn + u3*dt;

% Cartesian displacements
ua11(1,:,:) = real(SpHReconst(u1,Yt));
ua11(2,:,:) = real(SpHReconst(u2,Yt));
ua11(3,:,:) = real(SpHReconst(u3,Yt));


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! switching vel here to just stress
% ut = A2\bb2;

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

% Integral of f over surface
fintx = real(bb(1:3:end));
fintx = reshape(fintx,np,nt)';
finty = real(bb(2:3:end));
finty = reshape(finty,np,nt)';
fintz = real(bb(3:3:end));
fintz = reshape(fintz,np,nt)';

% q=quiver3(reshape(x1,[1,numel(x1)]),reshape(x2,[1,numel(x1)]),reshape(x3,[1,numel(x1)]),reshape(myf21,[1,numel(x1)]),reshape(myf22,[1,numel(x1)]),reshape(myf23,[1,numel(x1)]));

trcpnt(cts) = real(SpHReconst(x1mn,Ytrc));%x(1,floor(nt/2)+1,1);
trcvel(cts) = real(SpHReconst(u1,Ytrc));%ua(1,floor(nt/2)+1,1);
trcvelz(cts) = real(SpHReconst(u3,Ytrc));%ua(3,floor(nt/2)+1,1);

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

%vel plots
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
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if cts == 1 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append'); 
end 

% some output
disp([max(abs(x3mn-ix3mn)), real(u3(1)),max(max(max(abs(myf(:,:,:))))),cts]);

 end

















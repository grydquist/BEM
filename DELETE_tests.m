% Differences: get the deformed state info first. With that, calculate the
% velocity. Only need to get material constants once.

% NO VOLUME CORRECTION!

%% Material Constants

%notes: lam=1 corresponds to A = 4pi identity which falls out from the
%governing eq.

% Viscosity outside
mu = 1;
% Viscosity ratio
lam = 1;%/5;

% Non-dim numbers assume characteristic flow vel = 1 and radius =1

% Capillary number
Ca = 1/6;

% Dilatation ratio
Ed = 20;

% Bending ratio
Ebs = .015;

% Deformation resistance constants
B = 2/Ca;
C = Ed*B;
Eb = 2*B*Ebs;

% Spontaneous Curvature
c0 = 0;

% Total time steps
NT = 20;
dt = 0.005;

% Velocity and gradient
dU = [0,0,1;0,0,0;.0,0,0];
%% Reference Shape Definition

% Flag for later...
refflag = true;

% Order of the force and velocity transforms
p = 6;

% To de-alias the force, we need to get a finer grid. The factor is fali
fali = 4;

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
% Associated legendre polynomials at values of theta
mylegs = myleg(q,cos(thtf));
mylegsc= myleg(p,cos(tht));

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmTNew(p,th,ph);

% Harmonics evaluated at fine grid
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

x1mn = SpT(Yt,x1,th,ph);%0.5*SpT(Yt,x1,th,ph);
x2mn = SpT(Yt,x2,th,ph);%0.5*SpT(Yt,x2,th,ph);
x3mn = SpT(Yt,x3,th,ph);

% [x1mn,x2mn,x3mn] = BactCoeffs(Yt,th,ph);
x1mn(abs(x1mn)<1e-12) = 0;
x2mn(abs(x2mn)<1e-12) = 0;
x3mn(abs(x3mn)<1e-12) = 0;

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
ffs = myf;
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

gnR = zeros(2,2,ntf,npf);
dgtR = gnR;

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
Atot = zeros(3*nt*np, 3*ftot + 6);
At = zeros(3,3);
v = At;
b = zeros(3*nt*np,1);
bt = zeros(3,1);
thet = zeros(nt,np);
phit = thet;
A2 = zeros(3*ftot + 6,3*ftot + 6);
b2 = zeros(3*ftot,1);
Jg  = zeros(nt,np);
dxtg = zeros(3,1);
dxpg = dxtg;
xcg = zeros(3,nt,np);
ua = zeros(3,nt,np);
ua11 = ua;

%% Time stepping loop
cts = 1;
%% Calculation of residual force at convected points

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
        
%       Get the derivative in theta and add in
        dxt(1,:,:) = squeeze(dxt(1,:,:)) + f1*td1;
        dxt(2,:,:) = squeeze(dxt(2,:,:)) + f2*td1;
        dxt(3,:,:) = squeeze(dxt(3,:,:)) + f3*td1;
        
%       Get the derivative in phi and add in
        dxp(1,:,:) = squeeze(dxp(1,:,:)) + f1*1i*m*Ym;
        dxp(2,:,:) = squeeze(dxp(2,:,:)) + f2*1i*m*Ym;
        dxp(3,:,:) = squeeze(dxp(3,:,:)) + f3*1i*m*Ym;
    end
end
% These shouldn't be imaginary, but some creeps in        
dxt   = real(dxt);

dxp   = real(dxp);


% Geometry info
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

    end
end


%% Calculation of velocity constants via fluid problem

% Calculate the integral
ip = 0;
ic = 0;
A(:) = 0;
b(:) = 0;
nkg = nk;
At3 = zeros(3,3,3);
At2 = zeros(3,3);
AtG = At2;

% First loop: Inner integrals at Gauss points
for i = 1:nt
    for j = 1:np
        
%       Total Galerkin mode count
        ic = ic+1;
        
%       Velocity at colloc point
        Uc = dU*x(:,i,j);
        
%       Rotation Matrix
        t1 = [cos(phi(j)),-sin(phi(j)),0;sin(phi(j)),cos(phi(j)),0;0,0,1];    
        t2 = [cos(-tht(i)),0,sin(-tht(i));0,1,0;-sin(-tht(i)),0,cos(-tht(i))];
        t3 = [cos(-phi(j)),-sin(-phi(j)),0;sin(-phi(j)),cos(-phi(j)),0;0,0,1];

        Tx = t1*t2*t3;
          
%       Rotated harmonic coefficients of geometry
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
                dxtg(:) = 0;
                dxpg(:) = 0;
                
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
                        dxtg(1) = dxtg(1) + x1mnr(it)*td1;
                        dxtg(2) = dxtg(2) + x2mnr(it)*td1;
                        dxtg(3) = dxtg(3) + x3mnr(it)*td1;
                        
                        Ymn = m*Ypcur(im);
                        dxpg(1) = dxpg(1) + x1mnr(it)*Ymn;
                        dxpg(2) = dxpg(2) + x2mnr(it)*Ymn;
                        dxpg(3) = dxpg(3) + x3mnr(it)*Ymn;
                    end
                end
%               Now we can get the Jacobian here, and the normal vector
                
%               Inward normal
                nkg(:,i2,j2) = real(cross(dxtg,dxpg));
                nkg(:,i2,j2) = -nkg(:,i2,j2)./norm(nkg(:,i2,j2));
                
%               Jacobian (area element) via fundamental forms
                Jg(i2,j2) = real(sqrt(dot(dxtg,dxtg) ... 
                          * dot(dxpg,dxpg) ...
                          - dot(dxtg,dxpg)^2));
            end
        end
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        Ypt = SpHarmTNew(p,thet,phit);
        
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
                
                At3(:,:,:) = 0;
                AtG(:) = 0;
%               Calculate the integral part
                for ig = 1:nt
                    for jg = 1:np
                        r = (xcr-xcg(:,ig,jg));
                        T = newT(r);
%                         v = T(:,:,1)*nkg(1,ig,jg) + T(:,:,2)*nkg(2,ig,jg) + T(:,:,3)*nkg(3,ig,jg);
%                         At = At + v*Jg(ig,jg)*ws(ig)*Y(ig,jg)*dphi;
                        At3 = At3 + T*Jg(ig,jg)*ws(ig)*Y(ig,jg)*dphi;
                        v = Gij(r);
                        AtG = AtG + v*Y(ig,jg)*Jg(ig,jg)*ws(ig)/8/pi;
                    end
                end
               
%                 At = squeeze(At3(:,1,:)*nk(1,i,j) + At3(:,2,:)*nk(2,i,j) + At3(:,3,:)*nk(3,i,j));
%                 
%                 At = 3/4/pi*At + Yg(i,j)*eye(3);
%  
                col = 3*(ih-1)+1;
                A(row:row+2,col:col+2) = AtG;
            end
        end
        b(row:row+2) = Uc;
        
    end
end

% Levi-civita symbol
epsijk = zeros(3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            epsijk(i,j,k) = 0.5*(i-j)*(j-k)*(k-i);
        end
    end
end
% Integral constraint and full matrix
Atot(:,7:3*ftot + 6) = A;
ic  = 0;
for i = 1:nt
    for j = 1:np
        ic = ic+1;
        row = 3*(ic-1)+1;
%       Translational velocity
        Atot(row:row+2,1:3) = eye(3);
        
%       Rotational velocity
        xcr = x(:,i,j);
        mat = zeros(3,3);
        for k = 1:3
            mat = mat + squeeze(xcr(k)*epsijk(:,:,k));
        end
        Atot(row:row+2,4:6) = mat;
    end
end

%% Outer, Galerkin integral (done over a sphere!)
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
        
        Y = squeeze(Ypcur(im,:,:));
        
%       Loop over inner harmonics (a constant value here is a column)
        im2 = 0;
        for n2 = 0:p
            for m2 = -n2:n2
                im2 = im2+1;
                col = 3*im2 - 2 + 6;
%               Loop over Gauss points (these sum to get one row/col value)
                At(:) = 0;
                ii = 0;
                for i = 1:nt
                    for j = 1:np
                        ii = ii+1;
%                       Get value of inner integral at Gauss point
                        v = Atot(3*ii-2:3*ii,3*im2-2 + 6:3*im2 + 6);
%                       Integrate with that value
                        At = At + v*conj(Y(i,j))*wg(i)*dphi;
                    end
                end
                
                A2(row:row+2,col:col+2) = At;
            end
        end
        b2(row:row+2) = 0;
    end 
end

% Now for just the velocities
it = 0;
for n = 0:p
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        it = it+1;
        row = 3*it - 2;
        
        Y = squeeze(Ypcur(im,:,:));
        At(:) = 0;
        ii = 0;
        for i = 1:nt
            for j = 1:np
                ii = ii+1;
        %       Get value of inner integral at Gauss point
                v = Atot(3*ii-2:3*ii,1:3);
        %       Integrate with that value
                At = At + v*conj(Y(i,j))*wg(i)*dphi;
            end
        end
        
        A2(row:row+2,1:3) = At;
        
        At(:) = 0;
        ii = 0;
        for i = 1:nt
            for j = 1:np
                ii = ii+1;
        %       Get value of inner integral at Gauss point
                v = Atot(3*ii-2:3*ii,4:6);
        %       Integrate with that value
                At = At + v*conj(Y(i,j))*wg(i)*dphi;
            end
        end
        A2(row:row+2,4:6) = At;
    end
end

% Integral constraints, force then torque free
it = 0;
for n = 0:p
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        it = it+1;
        col = 3*it - 2;
        
        Y = squeeze(Ypcur(im,:,:));
        At(:) = 0;
        for i = 1:nt
            for j = 1:np
                At = At + Y(i,j)*wg(i)*dphi*eye(3);
            end
        end
        A2(end-5:end-3,col+6:col+2+6) = At;
        
    end
end

it = 0;
for n = 0:p
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        it = it+1;
        col = 3*it - 2;
        
        Y = squeeze(Ypcur(im,:,:));
        At(:) = 0;
        for i = 1:nt
            for j = 1:np
                xcr = x(:,i,j);
                mat = zeros(3,3);
                for k = 1:3
                    mat = mat + squeeze(xcr(k)*epsijk(:,:,k));
                end
                At = At + Y(i,j)*wg(i)*dphi*mat;
            end
        end
        A2(end-2:end,col+6:col+2+6) = At;
        
    end
end



% A2(end-2:end,:) = 0;
% A2(end-2,1:3:end) = 1;
% A2(end-1,2:3:end) = 1;
% A2(end-0,3:3:end) = 1;
b2(end-2:end) = 1;

ut = A2\b2;
ut(abs(ut)<1e-12) = 0;

u1 = ut(1:3:end);
u2 = ut(2:3:end);
u3 = ut(3:3:end);

f1 = SpHReconst(u1,Yt);
f2 = SpHReconst(u2,Yt);
f3 = SpHReconst(u3,Yt);







 function T = newT(r)
T = zeros(3,3,3);
rm = norm(r); 

 for i = 1:3
     for j = 1:3
         for k = 1:3
             T(i,j,k) = r(i)*r(j)*r(k)/rm^5;
         end
     end
 end
 
 
 end










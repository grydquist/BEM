load fmns_cvg.mat
x1mn = zeros(17*17,1);
x3mn = x1mn;
x1mn(2) = 1.437475235320810;
x1mn(4) = -1.437475235320810;
x2mn = x1mn*1i;
x2mn(4) = -x2mn(4);
x3mn(3) = 2.051708557987912;

% Max p is 16
p = 4;
nt =(p+1);
np = 2*nt;
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
[ph,th] = meshgrid(phi,tht);
thet = zeros(nt,np);
phit = thet;
Yt = SpHarmTNew(p,th,ph);
xcg = zeros(3,nt,np);
Nmyf = zeros(3,nt,np);
vG = zeros(3,3,nt,np);
Jg = zeros(nt,np);

ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end


Ytd1_c = zeros(nt,np,(p+1)^2);

it = 0;
im = 0;
for n = 0:p
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
        Ytd1_c(:,:,it)= ThetDer(Yt,ph,n,m,1);
    end
end

myfacs = zeros(2*(16)+1,1);
for i = 0:2*(16)
    myfacs(i+1) = factorial(i);
end

thint = 1;
phint = 1;
Ytcrint = SpHarmT(p,0,0,myfacs);

%       Rotation Matrix
t1 = [cos(phint),-sin(phint),0;sin(phint),cos(phint),0;0,0,1];
t2 = [cos(-thint),0,sin(-thint);0,1,0;-sin(-thint),0,cos(-thint)];
t3 = [cos(-phint),-sin(-phint),0;sin(-phint),cos(-phint),0;0,0,1]; 
Tx = t1*t2*t3;

%       Rotated harmonic coefficients of geometry

x1mnr = SpHRot(x1mn,phint,-thint,0-phint, myfacs);
x2mnr = SpHRot(x2mn,phint,-thint,0-phint, myfacs);
x3mnr = SpHRot(x3mn,phint,-thint,0-phint, myfacs);

%       Nonrotated locations of rotated grid points.        
xcg(1,:,:) = SpHReconst(x1mnr,Yt);
xcg(2,:,:) = SpHReconst(x2mnr,Yt);
xcg(3,:,:) = SpHReconst(x3mnr,Yt);
xcg = real(xcg);

%       The location of the new pole, should just be location of
%       point at current theta/phi
xcr(1) = SpHReconst(x1mnr,Ytcrint);
xcr(2) = SpHReconst(x2mnr,Ytcrint);
xcr(3) = SpHReconst(x3mnr,Ytcrint);
xcr = real(xcr);

%       Gauss points after rotation, in non-rotated reference.
for i2 = 1:nt
    for j2 = 1:np
%               Gauss point in rotated reference frame
        gp = [sin(tht(i2))*cos(phi(j2)); sin(tht(i2))*sin(phi(j2)); cos(tht(i2))];
        gp = Tx'*gp;
        phit(i2,j2) = atan2(gp(2),gp(1));
%               Because it can sometimtes be gp(3) - 1 = 2e-16
        if gp(3)> 1; gp(3) =  1;end
        if gp(3)<-1; gp(3) = -1;end
        thet(i2,j2) = acos(gp(3));

%               Some quick partial derivatives at locations in nonrotated
        dxtg = zeros(3,1);
        dxpg = zeros(3,1);

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
%       Now we can get the Jacobian here
%       Jacobian (area element) via fundamental forms
        Jg(i2,j2) = real(sqrt(dot(dxtg,dxtg) ... 
                  * dot(dxpg,dxpg) ...
                  - dot(dxtg,dxpg)^2));

%               Calculate kernels at rotated GPs
        r = xcr - xcg(:,i2,j2)';
        vG(:,:,i2,j2) = Gij(r);
    end
end
Ypt = SpHarmTNew(16,thet,phit);

%       Bookkeeping
bt(:) = 0;

%       Forces on rotated grid.
Nmyf(1,:,:) = real(SpHReconst(myfS1,Ypt));
Nmyf(2,:,:) = real(SpHReconst(myfS2,Ypt));
Nmyf(3,:,:) = real(SpHReconst(myfS3,Ypt));

for ig = 1:nt
    for jg = 1:np
        bt = bt + vG(:,:,ig,jg)*Nmyf(:,ig,jg)*Jg(ig,jg)*ws(ig)*dphi;
    end
end

%% Other integrals method

% Nonrotated locations of rotated grid points.        
xcg(1,:,:) = SpHReconst(x1mn,Yt);
xcg(2,:,:) = SpHReconst(x2mn,Yt);
xcg(3,:,:) = SpHReconst(x3mn,Yt);
xcg = real(xcg);

for i2 = 1:nt
    for j2 = 1:np

%       Some quick partial derivatives at locations in nonrotated
        dxtg = zeros(3,1);
        dxpg = zeros(3,1);

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
                dxtg(1) = dxtg(1) + x1mn(it)*td1;
                dxtg(2) = dxtg(2) + x2mn(it)*td1;
                dxtg(3) = dxtg(3) + x3mn(it)*td1;

                Ymn = m*Ypcur(im);
                dxpg(1) = dxpg(1) + x1mn(it)*Ymn;
                dxpg(2) = dxpg(2) + x2mn(it)*Ymn;
                dxpg(3) = dxpg(3) + x3mn(it)*Ymn;
            end
        end
%               Now we can get the Jacobian here

%               Jacobian (area element) via fundamental forms
        Jg(i2,j2) = real(sqrt(dot(dxtg,dxtg) ... 
                  * dot(dxpg,dxpg) ...
                  - dot(dxtg,dxpg)^2));

%               Calculate kernels at rotated GPs
        r = xcr - xcg(:,i2,j2)';
        vG(:,:,i2,j2) = Gij(r);
    end
end


Ypt = SpHarmTNew(16,th,ph);

%       Bookkeeping
bt2(:) = 0;

%       Forces on rotated grid.
Nmyf(1,:,:) = real(SpHReconst(myfS1,Ypt));
Nmyf(2,:,:) = real(SpHReconst(myfS2,Ypt));
Nmyf(3,:,:) = real(SpHReconst(myfS3,Ypt));

for ig = 1:nt
    for jg = 1:np
        bt2 = bt2 + vG(:,:,ig,jg)*Nmyf(:,ig,jg)*Jg(ig,jg)*wg(ig)*dphi;
    end
end

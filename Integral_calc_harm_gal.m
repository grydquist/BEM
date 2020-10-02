%% Trying to use BEM to calculate stress on objects in Stokes flow
% Most helpful here was probably Sorgentone:
% https://www.sciencedirect.com/science/article/pii/S0021999118300433
% And also Veerapani: boundary
mu = 1;
% U = [0;0;1];
U = [0;0;0];
nsd = 3;
dU = [0,0.5,0;0.5,0,0;0,0,0];
% dU = [0,.1,.1;.5,0,.1;.2,.3,0];
% dU = [0,0,0;0,0,0;0,0,0];
%% Actual code - Spherical harmonics creation

% Order of harmonics
p = 3;
% How many coefficients per dimension?
ftot = (p+1)^2;

% Our grid to evaluate integrals over the surface. We use a trapezoidal
% rule of equidistant points in phi, and Gauss points in theta
gf = 1;
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);
Yt = SpHarmT(p,th,ph);

% Weights for rotation 
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:gf*p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

%% Pre-processing required each time step

% Geometric coefficients
% xmn = [1/.2821,0.25,0,-.25];
% xmn = [1/0.2821,0,0,0];
% xmn = [1/0.2821,-.1 + .1i,.1,.1 + .1i,-.3,.2i,.2,.2i,-.3];
xmn = [4.5464,0,0,0,0,.1i,-0.6582,.1i,0];
% xmn = [4.5471,0,0,0,0,0,-0.6559,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.2132,0,0,0,0];
% xmn =1/0.2821;
% Order of geometric parameterization
pxmn = sqrt(length(xmn))-1;

% Harmonics at north pole
Ytcr= SpHarmT(pxmn,0,0);

% Points, normal vectors, Jacobian at Gauss points, not needed at this time
[x,nk,J] = SpHDer(xmn,Yt,th,ph);

%% Calculate the integral
ip = 0;
ic = 0;
A = zeros(3*nt*np, 3*ftot);
At = zeros(3,3);
v = At;
b = zeros(3*nt*np,1);
bt = zeros(3,1);
thet = zeros(nt,np);
phit = thet;

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
%               SpHarms order n, degreem eval'd at integ points
                Y = squeeze(Ypcur(im,:,:));
                At(:) = 0;
                for ig = 1:nt
                    for jg = 1:np
                        r = Tx'*(xcr-xcg(:,ig,jg));
                        v = Gij(r);
                        At = At + v*Y(ig,jg)*Jg(ig,jg)*ws(ig);
%                       Only need to calc B once per colloc point
                        if(n==0)
                            v = Tij(r,Tx'*nkg(:,ig,jg));
                            bt = bt + v*(U + dU*(Tx'*xcg(:,ig,jg)))*Jg(ig,jg)*ws(ig);
                        end
                    end
                end
 
                col = 3*(ih-1)+1;
                
%               Integral calc'd for colloc/harm combo. Put in A
                A(row:row+2,col:col+2) = A(row:row+2,col:col+2) + At*dphi;
            end
        end
%       Forcing integral calc'd for a colloc point! Put in B
        b(row:row+2) = b(row:row+2) + bt*dphi - 4*pi*Uc;
    end
end

%% Outer, Galerkin integral (done over a sphere!)
A2 = zeros(3*ftot,3*ftot);
b2 = zeros(3*ftot,1);
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
         
% AAc =A2;
% AAc(:,4:6) = real(A2(:,7:9));
% AAc(:,7:9) = real(A2(:,4:6)-A2(:,10:12));
% AAc(:,10:12) = imag(A2(:,4:6)+A2(:,10:12));
% AAc = real(AAc);
% AAc(abs(AAc)<.001) = 0;
% 
% rAc2 = zeros(12,24);
% iAc2 = rAc2;
% rAc2(:,1:2:end) = real(A2);
% rAc2(:,2:2:end) = -imag(A2);
% iAc2(:,1:2:end) = imag(A2);
% iAc2(:,2:2:end) = real(A2);
% rAc2(abs(rAc2)<.001) = 0;
% iAc2(abs(iAc2)<.001) = 0;
% 
% rAc2(:,2:2:6)=[];
% iAc2(:,2:2:6)=[];
% rAc2(:,10:2:14) = [];
% iAc2(:,10:2:14) = [];
% rAc2(:,4:2:8) = rAc2(:,4:2:8) - rAc2(:,13:2:17);
% iAc2(:,4:2:8) = iAc2(:,4:2:8) - iAc2(:,13:2:17);
% rAc2(:,5:2:9) = rAc2(:,5:2:9) + rAc2(:,14:2:18);
% iAc2(:,5:2:9) = iAc2(:,5:2:9) + iAc2(:,14:2:18);
% rAc2(:,13:18) = [];
% iAc2(:,13:18) = [];


[Us,S,V] = svd(A2);
Si = 1./S;
Si(Si>1e10)=0;
% Si(end) = Si(end-1,end-1);
Si(end) = 0;
Ai = V*Si*Us';
% f =A\b; 
f =Ai*b2;
% f = A2\b2;

%% Reconstruction
fa = zeros(nsd,nt,np);
pp = zeros(nt,np);
ind = 0;
for n = 0:p
%   All spherical harmonics of order n evaluated at integ points
    Ypcur = Yt{n+1};
    im = 0;
    for m = -n:n
        im = im+1;
        Y = Ypcur(im,:,:);
%         if(m~=0); Y2 = Ypcur(im-2*m,:,:); end
        for d = 1:nsd
            ind = ind+1;
%           SpHarms order n, degree m eval'd at integ points
            fa(d,:,:) = fa(d,:,:) + Y*f(ind);
%             if(m~=0)
%                 fa(d,:,:) = fa(d,:,:) + Y2*((-1)^m*conj(f(ind)));
%             end
        end
    end
end
fa = real(fa);
for i =1:nt
    for j =1:np
        pp(i,j) = real(fa(1,i,j)*nk(1,i,j) + fa(2,i,j)*nk(2,i,j) + fa(3,i,j)*nk(3,i,j));
    end
end


%% Plots
% plot(tht,pp(:,1));
% 
% 
% x1 = zeros(1,ntc*npc);
% y1 = x1;
% z1 = x1;
% cnt= 0;
% for i = 1:ntc
%     for j = 1:npc
%         cnt = cnt+1;
%         x1(cnt) = xc(1,i,j);
%         y1(cnt) = xc(2,i,j);
%         z1(cnt) = xc(3,i,j);
%     end
% end
% 
% 
x2 = zeros(1,nt*np);
y2 = x2;
z2 = x2;
pf = x2;
cnt= 0;
for i = 1:nt
    for j = 1:np
        cnt = cnt+1;
        x2(cnt) = real(x(1,i,j));
        y2(cnt) = real(x(2,i,j));
        z2(cnt) = real(x(3,i,j));
        pf(cnt) = pp(i,j);
    end
end
figure;
scatter3(x2,y2,z2,100,pf,'filled')
pbaspect([1,1,1])
axis([-1,1,-1,1,-1,1])
xlabel('X')
ylabel('Y')
% 
% % Ynm = SpHarm(3,th,ph); 
% % Ymn = squeeze(Ynm(2,:,:));
% % 
% % [Xm,Ym,Zm] = sph2cart(ph, th-pi/2, real(Ymn).^2);
% % surf(Xm,Ym,Zm,'edgecolor','none')
% % 
% % Ymn =squeeze(Ypcur(1,:,:));
% % [Xm,Ym,Zm] = sph2cart(phit, thet-pi/2, real(Ymn).^2);
% % surf(Xm,Ym,Zm,'edgecolor','none')
% % scatter3(reshape(Xm,[1,numel(Xm)]),reshape(Ym,[1,numel(Xm)]),reshape(Zm,[1,numel(Xm)]))
% 
% 
% 
% % [Xm,Ym,Zm] = sph2cart(phit, thet-pi/2, ones(size(thet)));
% % scatter3(reshape(Xm,[1,numel(Xm)]),reshape(Ym,[1,numel(Xm)]),reshape(Zm,[1,numel(Xm)]),100, reshape(squeeze(Gt(:,:,i,j,1,2)),[1,numel(Zm)]), 'filled')
hold on
% Now for the velocities at the integration points
Ug = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        Ug(:,i,j) = U + dU*x(:,i,j);
    end
end
x = real(x);
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(Ug(1,:,:),[1,numel(x(1,:,:))]),reshape(Ug(2,:,:),[1,numel(x(1,:,:))]),reshape(Ug(3,:,:),[1,numel(x(1,:,:))]))
pbaspect([1,1,1])
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
% 
figure;
Yy = SpHarmT(pxmn,th,ph);
f1r = SpHReconst(xmn,Yy);
[Xm,Ym,Zm] = sph2cart(ph, th-pi/2, f1r);
Xm = real(Xm);
Ym = real(Ym);
Zm = real(Zm);
surf(Xm,Ym,Zm,pp,'edgecolor','none')
pbaspect([1,1,1])
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
%% Integral test
% Ag = zeros(3,3*ftot);
% ih=0;
% bt=0;
% Bt = zeros(3,3);
% tx = zeros(1,nt);
% ic = 0;
% % Loop over harmonics
% thet = squeeze(alla(1,1,1,:,:));
% phit = squeeze(alla(2,1,1,:,:));
% for n = 0:p
% %   All spherical harmonics of order n evaluated at integ points
% %     Ypcur = SpHarm(n,thet,phit);
%     Ypcur = SpHarm(n,th,ph);
%     im = 0;
%     for m = -n:n
%         ih = ih+1;
%         im = im+1;
% %       SpHarms order n, degreem eval'd at integ points
%         Y = squeeze(Ypcur(im,:,:));
%         At(:) = 0;
%         At2(:) = 0;
%         for ig = 1:np
%             for jg = 1:nt
% %                 r = squeeze(e(1,1,:,:))\([0,0,1]'-x(:,jg,ig));
%                 r = ([0,0,1]'-x(:,jg,ig));
%                 v = Gij(r);
%                 At = At + v*Y(jg,ig)*ws(jg)*dphi*Jg(jg,ig);
%                 Atr = real(At);
%                 Ati = imag(At);
% 
% %               Only need to calc B once per colloc point
%                 if(n==0)
%                     v = Tij(r,-x(:,jg,ig));
%                     bt = bt + v*Ug(:,jg,ig)*ws(jg)*dphi;
%                     
%                     Bt = Bt + v*ws(jg)*dphi;
%                     if(ig==1); tx(jg) = v(3,3);end
%                 end
%             end
%         end
%         col = 3*(ih-1)+1;
%         Ag(:,col:col+2) =(Atr + 1i*Ati);
%     end
% end
% bg = bt - 4*pi*[0;0;0];
% 
% % f2 = zeros(3*ftot,1);
% % im1 = 0;
% % im2 = 0;
% % for n = 0:p
% %     for m = -n:n
% %         if(m>=0);im1 = im1+1;end
% %         im2 = im2+1;
% %         if(m>=0)
% %             f2(3*im2 -2:3*im2) = f(3*im1-2:3*im1);
% %         else
% %             f2(3*im2 -2:3*im2) = (-1)^(m)*conj(f(3*(im1-2*m)-2:3*(im1-2*m)));
% %         end
% %     end
% % end
% 
% % gnm = Ag(3,3:3:end);
% % gnm = real(gnm);
% % gnm([abs(gnm)<0.0001])=[];
% % r = [0,0,1]'-x(:,1,1);
% % comp = Gij(r);
% % c2 = 0;
% % for n = 0:length(gnm)-1
% %     Yy = SpHarm(n,tht(1),phi(1));
% %     c2 = c2 + gnm(n+1)*Yy(n+1);
% % end
% L=zeros(3,1);
% for i = 1:nt
%     for j = 1:np
%         L = L + cross(fa(:,i,j),nk(:,i,j))*J(i,j)*dphi*wg(i);
%     end
% end
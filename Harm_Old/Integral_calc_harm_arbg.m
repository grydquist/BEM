%% Trying to use BEM to calculate stress on objects in Stokes flow
% Most helpful here was probably Sorgentone:
% https://www.sciencedirect.com/science/article/pii/S0021999118300433
% And also Veerapani: boundary
mu = 1;
% U = [0;0;1];
U = [0;0;0];
nsd = 3;
% dU = [0,0,0;0,0,0;0,0,0];
% dU = [0,.1,.1;.5,0,.1;.2,.3,0];
dU = [0,0.5,0;0.5,0,0;0,0,0];
%% Actual code - Spherical harmonics creation

% Order of harmonics
p = 4;
% How many coefficients per dimension?
ftot = (p+1)^2;
% Colloc points. Really bad way to do it, but just trying to get it working
ntc = p+1;
npc = p+1;

thtc = (0:ntc+1)/(ntc+1)*pi;
thtc = thtc(2:end-1)';
phic = (2*pi/npc/2):(2*pi/npc):(2*pi/npc)*npc-(2*pi/npc/2);
phic = phic-pi/4;

% Our grid to evaluate integrals over the surface. We use a trapezoidal
% rule of equidistant points in phi, and Gauss points in theta
gf = 1;
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1); % This seems to maybe take a lon time ??????
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);
[phc,thc] = meshgrid(phic,thtc);
Yt = SpHarmT(p,th,ph);

% Weird weights 
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:gf*p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

%% Pre-processing required each time step

% Geometric coefficients
% xmn = [1,0.25,0,-.25];
% xmn = [1/0.2821,0,0,0];
% xmn = [1/0.2821,-.1 + .1i,.1,.1 + .1i,-.3,.2i,.2,.2i,-.3];
% xmn = [4.5464,0,0,0,0,.1i,-0.6582,.1i,0];
% xmn = [4.5464,0,0,0,0,0,-0.6582,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.2178,0,0,0,0];
xmn =1/0.2821;
% Order of geometric parameterization
pxmn = sqrt(length(xmn))-1;

Ytc = SpHarmT(pxmn,thc,phc);
% Harmonics at north pole
Ytcr= SpHarmT(pxmn,0,0);

% We will also need the Jacobian over our GEOMETRY for each point.
% This is easy for our sphere test case, but in the future it should be the
% ||dr/dth X drd/dph||, where r = (x, y, z), and x =
% cos(thet)sin(thet) sum^(nm)anmYnm, etc. Likely will need to make good use
% of product rule. Also other formulae, see stack overflow/Veerapani

[x,nk,J] = SpHDer(xmn,Yt,th,ph);
xc = SpHDer(xmn,Ytc,thtc,phic); % Probably don't need to do the full Der!!!

%% Calculate the integral
ip = 0;
it = 0;
ic = 0;
A = zeros(3*ftot, 3*ftot);
At = zeros(3,3);
v = At;
b = zeros(3*ftot,1);
bt = zeros(3,1);
thet = zeros(nt,np);
phit = thet;

for i = 1:ntc
    for j = 1:npc
        
%       Collocation count
        ic = ic+1;
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        
%       Velocity at colloc point
        Uc = U + dU*xc(:,i,j);        
        
%       Rotation Matrix
        t1 = [cos(phic(j)),-sin(phic(j)),0;sin(phic(j)),cos(phic(j)),0;0,0,1];
        t2 = [cos(-thtc(i)),0,sin(-thtc(i));0,1,0;-sin(-thtc(i)),0,cos(-thtc(i))];
        t3 = [cos(-phic(j)),-sin(-phic(j)),0;sin(-phic(j)),cos(-phic(j)),0;0,0,1]; 
        Tx = t1*t2*t3;

%       Gauss points after rotation, in non-rotated reference. !!! likely a better way to do this
%       I'm not actually sure if this is event he right way to do this!!!!!
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
                thet(i2,j2) = acos(gp(3));
            end
        end
        
%       Rotated harmonic coefficients of geometry
%       Perhaps we don't need to do this? The Jacobian should be rotation
%       invariant and we rotate the normal vector back anyways?????????????
        xmnr = SpHRot(xmn,phic(j),-thtc(i),-phic(j));
        
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

% AAc =A;
% AAc(:,4:6) = A(:,7:9);
% AAc(:,7:9) = A(:,4:6)-A(:,10:12);
% AAc(:,10:12) = imag(A(:,4:6)+A(:,10:12));
[Us,S,V] = svd(A);
Si = 1./S;
Si(Si>1e10) = 0;
% Si((end-3*(p-1)):end,(end-3*(p-1)):end) = 0;
Ai = V*Si*Us';
% f =A\b; % LOOK OUT!!!!!!!!!!!! This still doesn't work... for radial sym?
f =Ai*b;

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


x1 = zeros(1,ntc*npc);
y1 = x1;
z1 = x1;
cnt= 0;
for i = 1:ntc
    for j = 1:npc
        cnt = cnt+1;
        x1(cnt) = xc(1,i,j);
        y1(cnt) = xc(2,i,j);
        z1(cnt) = xc(3,i,j);
    end
end


x2 = zeros(1,ntc*npc);
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

% Ynm = SpHarm(3,th,ph); 
% Ymn = squeeze(Ynm(2,:,:));
% 
% [Xm,Ym,Zm] = sph2cart(ph, th-pi/2, real(Ymn).^2);
% surf(Xm,Ym,Zm,'edgecolor','none')
% 
% Ymn =squeeze(Ypcur(1,:,:));
% [Xm,Ym,Zm] = sph2cart(phit, thet-pi/2, real(Ymn).^2);
% surf(Xm,Ym,Zm,'edgecolor','none')
% scatter3(reshape(Xm,[1,numel(Xm)]),reshape(Ym,[1,numel(Xm)]),reshape(Zm,[1,numel(Xm)]))



% [Xm,Ym,Zm] = sph2cart(phit, thet-pi/2, ones(size(thet)));
% scatter3(reshape(Xm,[1,numel(Xm)]),reshape(Ym,[1,numel(Xm)]),reshape(Zm,[1,numel(Xm)]),100, reshape(squeeze(Gt(:,:,i,j,1,2)),[1,numel(Zm)]), 'filled')
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
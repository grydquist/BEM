%% Trying to use BEM to calculate stress on objects in Stokes flow
% Most helpful here was probably Sorgentone:
% https://www.sciencedirect.com/science/article/pii/S0021999118300433
% And also Veerapani: boundary
mu = 1;
U = [0;0;0];
nsd = 3;
dU = [0,0,0;0,0,0;0,0,0];
dU = [0,.5,0;-.5,0,0;0,0,0];
s = rng;
%% Actual code - Spherical harmonics creation

% Order of harmonics !!! Gets it exactly right for even p??? Think I'm
% maybe somehow doing something wrong with odds
p = 2;
ff = (p+1)*(p+1);
% How many coefficients per dimension?
ftot = (p+1)^2;
% Colloc points
ntc = p+1;
npc = p+1;


thtc = (0:ntc+1)/(ntc+1)*pi;
thtc = thtc(2:end-1);
phic = (2*pi/npc/2):(2*pi/npc):(2*pi/npc)*npc-(2*pi/npc/2);
phic = phic-pi/4;

% Our grid to evaluate integrals over the surface. We use a trapezoidal
% rule of equidistant points in phi, and Gauss points in theta
gf = 2;
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);
Yt = cell(p+1,1);
for i = 0:p
    Yt{i+1} = SpHarm(i,th,ph);
end

% Weird weights 
ws = zeros(1,length(wg));
for i = 1:nt
    for j = 0:gf*p
        ws(i) = ws(i) + legendreP(j,xs(i))/cos(tht(i)/2);
    end
    ws(i) = ws(i)*wg(i);
end

%%%%% Also be sure to calculate derivatives of the individual spherical
%%%%% harmonics at each point. Will make life easier when calculating
%%%%% Jacobian and normal vector

%% Pre-processing required each time step

% We will also need the Jacobian over our GEOMETRY for each point.
% This is easy for our sphere test case, but in the future it should be the
% ||dr/dth X drd/dph||, where r = (x, y, z), and x =
% cos(thet)sin(thet) sum^(nm)anmYnm, etc. Likely will need to make good use
% of product rule
J = zeros(nt,np);
nk = zeros(nsd,nt,np);
% x = zeros(3,nt,np);
for i = 1:nt
    for j =1:np
        J(i,j) = sin(tht(i));
        nk(1,i,j) = sin(tht(i))^2*cos(phi(j));
        nk(2,i,j) = sin(tht(i))^2*sin(phi(j));
        nk(3,i,j) = sin(tht(i))*cos(tht(i));
        nk(:,i,j) = nk(:,i,j)/sqrt(sum(nk(:,i,j).^2));
    end
end
x = nk;
nk = -nk;

% Now let's find the collocation points/vel at colloc
xc = zeros(3,ntc,npc);
Uc = zeros(3,ntc,npc);
nc = zeros(3,ntc,npc);
e = zeros(ntc,npc,3,3);
% rng(s);
for i = 1:ntc
    for j =1:npc
%         thtc(i) = rand(1)*pi;
%         phic(j) = rand(1)*2*pi;
        
        xc(1,i,j) = sin(thtc(i))^2*cos(phic(j));
        xc(2,i,j) = sin(thtc(i))^2*sin(phic(j));
        xc(3,i,j) = sin(thtc(i))*cos(thtc(i));
        xc(:,i,j) = xc(:,i,j)/sqrt(sum(xc(:,i,j).^2));
%         if(i==1 && j==1)
%             xc(1,i,j) = sin(.1*pi)^2*cos(.1*pi);
%             xc(2,i,j) = sin(.1*pi)^2*sin(.1*pi);% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             xc(3,i,j) = sin(.1*pi)*cos(.1*pi);
%             xc(:,i,j) = xc(:,i,j)/sqrt(sum(xc(:,i,j).^2));
%         end
%         if(i==2 && j==2)
%             xc(1,i,j) = sin(1.2*pi)^2*cos(1.2*pi);
%             xc(2,i,j) = sin(1.2*pi)^2*sin(1.2*pi);% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             xc(3,i,j) = sin(1.2*pi)*cos(1.2*pi);
%             xc(:,i,j) = xc(:,i,j)/sqrt(sum(xc(:,i,j).^2));
%         end
%         if(i==2 && j==1)
%             xc(1,i,j) = sin(0.4*pi)^2*cos(1.8*pi);
%             xc(2,i,j) = sin(.4*pi)^2*sin(1.8*pi);% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             xc(3,i,j) = sin(.4*pi)*cos(1.8*pi);
%             xc(:,i,j) = xc(:,i,j)/sqrt(sum(xc(:,i,j).^2));
%         end
        Uc(:,i,j) = U + dU*xc(:,i,j);
%       Rotation Matrix
        t1 = [cos(phic(j)),-sin(phic(j)),0;sin(phic(j)),cos(phic(j)),0;0,0,1];
        t2 = [cos(-thtc(i)),0,sin(-thtc(i));0,1,0;-sin(-thtc(i)),0,cos(-thtc(i))];
        t3 = [cos(-phic(j)),-sin(-phic(j)),0;sin(-phic(j)),cos(-phic(j)),0;0,0,1];
        
%         if(i==1 && j==1)
%             t1 = [cos(.1*pi),-sin(.1*pi),0;sin(.1*pi),cos(.1*pi),0;0,0,1];
%             t2 = [cos(-.1*pi),0,sin(-.1*pi);0,1,0;-sin(-.1*pi),0,cos(-.1*pi)];% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             t3 = [cos(-.1*pi),-sin(-.1*pi),0;sin(-.1*pi),cos(-.1*pi),0;0,0,1];
%         end
%         if(i==2 && j==2)
%             t1 = [cos(1.2*pi),-sin(1.2*pi),0;sin(1.2*pi),cos(1.2*pi),0;0,0,1];
%             t2 = [cos(-1.2*pi),0,sin(-1.2*pi);0,1,0;-sin(-1.2*pi),0,cos(-1.2*pi)];% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             t3 = [cos(-1.2*pi),-sin(-1.2*pi),0;sin(-1.2*pi),cos(-1.2*pi),0;0,0,1];
%         end
%         if(i==1 && j==2)
%             t1 = [cos(1.8*pi),-sin(1.8*pi),0;sin(1.8*pi),cos(1.8*pi),0;0,0,1];
%             t2 = [cos(-.4*pi),0,sin(-.4*pi);0,1,0;-sin(-.4*pi),0,cos(-.4*pi)];% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             t3 = [cos(-1.8*pi),-sin(-1.8*pi),0;sin(-1.8*pi),cos(-1.8*pi),0;0,0,1];
%         end
        Tx = t1*t2*t3;
        e(i,j,:,:) = Tx;
    end
end

% Now for the velocities at the integration points
Ug = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        Ug(:,i,j) = U + dU*x(:,i,j);
    end
end


%% Calculate angles of Gauss points in rotated reference frame.
% If this works, this could be accomplished by instead just rotating the
% actual spherical harmonics, which would be cheaper and less memory
% intense
Ut = zeros(3,nt,np,ntc,npc);
alla = zeros(2,ntc,npc,nt,np);
gpt = zeros(1,2);
Jt = zeros(ntc,npc,nt,np);
for i1 = 1:ntc
    for j1 = 1:npc
%       We're at a collocation point. Now let's see what the angles are
%       pre-rotation
        for i2 = 1:nt
            for j2 = 1:np
%               Normal vector to current colloc point
                a = xc(:,i1,j1);
                Tx = squeeze(e(i1,j1,:,:));
%               Gauss point in rotated reference frame
                gp = [sin(tht(i2))*cos(phi(j2)), sin(tht(i2))*sin(phi(j2)), cos(tht(i2))];
                gp = gp*Tx;
%               Theta
                gpt(2) = atan2(gp(2),gp(1));
                gpt(1) = acos(gp(3));
                alla(:,i1,j1,i2,j2) = gpt;
%               Jacobian's gonna need some real work
%               Better way would probably just be to have the area element
%               in the rotated reference frame, as the Jacobian for sphere
%               to sphere is technically 1.

%               In time, let's change it to the area element, so we don't
%               even need to go through the sphere.
                Jt(i1,j1,i2,j2) = sin(tht(i2));
                Ut(:,i2,j2,i1,j1) = U + dU*gp';
            end
        end
    end
end

%% Calculate the integral
ip = 0;
it = 0;
ic = 0;
% Split these up into real and imaginary so I can do conjugate stuff
A = zeros(3*ff, 3*ff);
At = zeros(3,3);
v = At;
b = zeros(3*ff,1);
bt = zeros(3,1);


for i = 1:ntc
    for j = 1:npc
%       Collocation count
        ic = ic+1;
        
%       Bookkeeping
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
        
%       Points to recalculate spherical harmonics at
        thet = squeeze(alla(1,i,j,:,:));
        phit = squeeze(alla(2,i,j,:,:));
        
%       Velocity at colloc point
        Uc = U + dU*xc(:,i,j);
        
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
                        r = (squeeze(e(i,j,:,:)))\([0,0,1]'-x(:,ig,jg));
                        v = Gij(r); 
                        At = At + v*Y(ig,jg)*Jt(i,j,ig,jg)*ws(ig)*dphi;
%                       Only need to calc B once per colloc point
                        if(n==0)
                            v = Tij(r,-(squeeze(e(i,j,:,:)))\x(:,ig,jg));
                            bt = bt + v*(U + dU*(squeeze(e(i,j,:,:))\x(:,ig,jg)))*Jt(i,j,ig,jg)*ws(ig);
                        end
                    end
                end
 
                col = 3*(ih-1)+1;
                
%               Integral calc'd for colloc/harm combo. Put in A
                A(row:row+2,col:col+2) = A(row:row+2,col:col+2) + At;
            end
        end
%       Forcing integral calc'd for a colloc point! Put in B
        b(row:row+2) = b(row:row+2) + bt*dphi - 4*pi*Uc;
    end
end


AAc =A;
AAc(:,4:6) = A(:,7:9);
AAc(:,7:9) = A(:,4:6)-A(:,10:12);
AAc(:,10:12) = imag(A(:,4:6)+A(:,10:12));
[U,S,V] = svd(A);
Si = 1./S;
Si(Si>1e10) = 0;
% Si(end-2:end,end-2:end) = 0;
Ai = V*Si*U';
f =Ai*b;
% f =A\b;

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
        pp(i,j) = fa(1,i,j)*nk(1,i,j) + fa(2,i,j)*nk(2,i,j) + fa(3,i,j)*nk(3,i,j);
    end
end


%% Plots
plot(tht,pp(:,1));


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
        x2(cnt) = x(1,i,j);
        y2(cnt) = x(2,i,j);
        z2(cnt) = x(3,i,j);
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
quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(Ug(1,:,:),[1,numel(x(1,:,:))]),reshape(Ug(2,:,:),[1,numel(x(1,:,:))]),reshape(Ug(3,:,:),[1,numel(x(1,:,:))]))

%% Integral test
Ag = zeros(3,3*ftot);
ih=0;
bt=0;
Bt = zeros(3,3);
tx = zeros(1,nt);
ic = 0;
% Loop over harmonics
thet = squeeze(alla(1,1,1,:,:));
phit = squeeze(alla(2,1,1,:,:));
for n = 0:p
%   All spherical harmonics of order n evaluated at integ points
%     Ypcur = SpHarm(n,thet,phit);
    Ypcur = SpHarm(n,th,ph+3*pi/2);
    im = 0;
    for m = -n:n
        ih = ih+1;
        im = im+1;
%       SpHarms order n, degreem eval'd at integ points
        Y = squeeze(Ypcur(im,:,:));
        At(:) = 0;
        At2(:) = 0;
        for ig = 1:np
            for jg = 1:nt
%                 r = squeeze(e(1,1,:,:))\([0,0,1]'-x(:,jg,ig));
                r = ([0,0,1]'-x(:,jg,ig));
                v = Gij(r);
                At = At + v*Y(jg,ig)*ws(jg)*dphi;
                Atr = real(At);
                Ati = imag(At);

%               Only need to calc B once per colloc point
                if(n==0)
                    v = Tij(r,-x(:,jg,ig));
                    bt = bt + v*Ug(:,jg,ig)*ws(jg)*dphi;
                    
                    Bt = Bt + v*ws(jg)*dphi;
                    if(ig==1); tx(jg) = v(3,3);end
                end
            end
        end
        col = 3*(ih-1)+1;
        Ag(:,col:col+2) =(Atr + 1i*Ati);
    end
end
bg = bt - 4*pi*[0;0;0];

% f2 = zeros(3*ftot,1);
% im1 = 0;
% im2 = 0;
% for n = 0:p
%     for m = -n:n
%         if(m>=0);im1 = im1+1;end
%         im2 = im2+1;
%         if(m>=0)
%             f2(3*im2 -2:3*im2) = f(3*im1-2:3*im1);
%         else
%             f2(3*im2 -2:3*im2) = (-1)^(m)*conj(f(3*(im1-2*m)-2:3*(im1-2*m)));
%         end
%     end
% end

% gnm = Ag(3,3:3:end);
% gnm = real(gnm);
% gnm([abs(gnm)<0.0001])=[];
% r = [0,0,1]'-x(:,1,1);
% comp = Gij(r);
% c2 = 0;
% for n = 0:length(gnm)-1
%     Yy = SpHarm(n,tht(1),phi(1));
%     c2 = c2 + gnm(n+1)*Yy(n+1);
% end
L=zeros(3,1);
for i = 1:nt
    for j = 1:np
        L = L + cross(fa(:,i,j),nk(:,i,j))*J(i,j)*dphi*wg(i);
    end
end
% Most helpful here was probably Sorgentone:
% https://www.sciencedirect.com/science/article/pii/S0021999118300433
% And also Veerapani: boundary
mu = 1;
U = [0;0;1];
nsd = 3;
dU = [0,0,0;0,0,0;0,0,0];
%% Actual code - Spherical harmonics creation

% Order of harmonics
p = 1;

ff = (p+1)*(p+2)/2;
% How many coefficients per dimension?
ftot = (p+1)^2;
% Colloc points. Really bad way to do it, but just trying to get it working
ntc = p+1;
npc = p+1;


thtc = (0:ntc+1)/(ntc+1)*pi;
thtc = thtc(2:end-1);
phic = (0:npc)/(npc)*2*pi;
phic = phic(1:end-1);

% Our grid to evaluate integrals over the surface. We use a trapezoidal
% rule of equidistant points in phi, and Gauss points in theta
np = 16*(p+2);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = 8*(p+2);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);
Yt = cell(p+1,1);
for i = 0:p
    Yt{i+1} = SpHarm(i,th,ph);
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
x = zeros(3,nt,np);
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

% Now let's find the collocation points
xc = zeros(3,ntc,npc);
Uc = zeros(3,ntc,npc);
for i = 1:ntc
    for j =1:npc
        xc(1,i,j) = sin(thtc(i))^2*cos(phic(j));
        xc(2,i,j) = sin(thtc(i))^2*sin(phic(j));
        xc(3,i,j) = sin(thtc(i))*cos(thtc(i));
        xc(:,i,j) = xc(:,i,j)/sqrt(sum(xc(:,i,j).^2));
        Uc(:,i,j) = U + dU*xc(:,i,j);
    end
end

% Now for the velocities at the integration points
Ug = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        Ug(:,i,j) = U + dU*x(:,i,j);
    end
end

%% Calculate the integral
a=0;

ip = 0;
it = 0;
ic = 0;
% Split these up into real and imaginary so I can do conjugate stuff
A = zeros(3*ff, 3*ff);
Ar = A;
Ai = A;
Atr = zeros(3,3);
Ati = Atr;
At = Atr;
At2 = At;
v = Atr;
b = zeros(3*ff,1);
bt = zeros(3,1);

for i = 1:npc
    for j = 1:ntc
%       Collocation count
        ic = ic+1;
%       Hacky way to have colloc points equal number of unknowns
        if(ic>ff); break; end
        row = 3*(ic-1)+1;
        bt(:) = 0;
        ih = 0;
%       Loop over harmonics
        for n = 0:p
%           All spherical harmonics of order n evaluated at integ points
            Ypcur = Yt{n+1};
            im = n; %???
            for m = 0:n
                ih = ih+1;
                im = im+1;
%               SpHarms order n, degreem eval'd at integ points
                Y = squeeze(Ypcur(im,:,:));
%               Corresponding negative harmonic
                if(m~=0); Y2 = squeeze(Ypcur(im-2*m,:,:)); end
                At(:) = 0;
                At2(:) = 0;
                for ig = 1:np
                    for jg = 1:nt
                        r = xc(:,j,i)-x(:,jg,ig);
                        v = Gij(r);
                        At = At + v*Y(jg,ig) ...
                            *J(jg,ig)*wg(jg)*dphi/sin(tht(jg));
                        Atr = real(At);
                        Ati = imag(At);
%                       Only need to calc B once per colloc point
                        if(n==0)
                            bt = bt + Tij(r,nk(:,jg,ig))*Ug(:,jg,ig) ... 
                                *J(jg,ig)*wg(jg)*dphi/sin(tht(jg));
                        end
                        
%                       If m isn't 0, add the corresponding negative f in 
                        if(m~=0)
                            At2 = At2 + v*Y2(jg,ig) ...
                            *J(jg,ig)*wg(jg)*dphi/sin(tht(jg));
                            Atr = Atr + (-1)^m*real(At2);
                            Ati = Ati - (-1)^m*imag(At2);
                        end
                    end
                end
                col = 3*(ih-1)+1;
                
%               Integral calc'd for colloc/harm combo. Put in A
                Ar(row:row+2,col:col+2) = Ar(row:row+2,col:col+2) + Atr;
                Ai(row:row+2,col:col+2) = Ai(row:row+2,col:col+2) + Ati;
            end
        end
%       Forcing integral calc'd for a colloc point! Put in B
        b(row:row+2) = b(row:row+2) + bt - 4*pi*Uc(:,j,i);
    end
end
A = Ar + 1i*Ai;
f =A\b;

%% Reconstruction
fa = zeros(nsd,nt,np);
pp = zeros(nt,np);
ind = 0;
for n = 0:p
%   All spherical harmonics of order n evaluated at integ points
    Ypcur = Yt{n+1};
    im = n;
    for m = 0:n
        im = im+1;
        Y = Ypcur(im,:,:);
        if(m~=0); Y2 = Ypcur(im-2*m,:,:); end
        for d = 1:nsd
            ind = ind+1;
%           SpHarms order n, degree m eval'd at integ points
            fa(d,:,:) = fa(d,:,:) + Y*f(ind);
            if(m~=0)
                fa(d,:,:) = fa(d,:,:) + Y2*((-1)^m*conj(f(ind)));
            end
        end
    end
end
%fa = real(fa);
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

% Ynm = SpHarm(3,th,ph); 
% Ymn = squeeze(Ynm(2,:,:));
% 
% [Xm,Ym,Zm] = sph2cart(ph, th-pi/2, real(Ymn).^2);
% surf(Xm,Ym,Zm,'edgecolor','none')



%% Integral test
% Yy = Yt{3};
% Y1 = squeeze(Yy(5,:,:));
% Yy = Yt{1};
% Y2 = squeeze(Yy(1,:,:));
% I = 0;
% 
% for i2 = 1:np
%     for j2 = 1:nt
% for ig = 1:np
%     for jg = 1:nt
%         I = I + Y1(jg,ig)*conj(Y1(j2,i2))...
%                             *J(jg,ig)*wg(jg)*dphi/sin(tht(jg))...
%                             *J(j2,i2)*wg(j2)*dphi/sin(tht(j2));
%     end
% end
%     end
% end


% Calculate Gij everywhere first!

% % First let's calculate the inner integral at each Gauss point
% for n = 0:p-1
%     for m = -n:n
%         for l=n:p
%             for o=-l:l
%                 for i=1:np-1
%                     for j = 1:nt-1
%                         for k = i:np
%                             for z=j:nt
%                                 a = a+1;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end


% Outer, Galerkin integral
% Loop over integration points
% for ip1 = 1:np
%     for jt2 = 1:nt
% %       Loop over outer Galerkin Spherical harmonics
%         for n1 = 0:p
%             for m1 = -n1:n1
% %               Inner, basis function integrals
% %               Loop over integration points
%                 for ip2 = 1:np
%                     for jt2 = 1:nt
% %                       Loop over inner basis spherical harmonics
%                         for n2 = 0:p
%                             for m2 = -n2:n2
%                                 a = a+1;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end



%%
% % Plot tests

% % Some GOOD plots of pressures and whatnot at Gauss points


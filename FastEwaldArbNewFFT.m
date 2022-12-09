%% For arbitrary/random points
%%% See Greengard and Lee (via tornberg) for speeding up

% First, define the lattice (columns vectors)
bv = [1 0 0; 0 1 0; 0 0 1]*1;
% bv=[1.0, 0.31, 5E-002
%     0.1,  1.0,  0.2
%     0.12,0.19,  1.0]';
bvx = bv(:,1);
bvy = bv(:,2);
bvz = bv(:,3);
x = [.5;.75,;.5];

% reciprocal basis vectors
kv = bv;
tau = dot(cross(bvx,bvy),bvz);
kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3));
kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1));
kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2));
kvx = kv(:,1);
kvy = kv(:,2);
kvz = kv(:,3);

bxs = 1;


% % Random x's
% pts1 = 10;
% pts = pts1^2;
% x0 = zeros(3,pts);
% f = zeros(3,pts);
% rng(69420);
% for i = 1:pts
% %   forces alternate up and down
%     f(1,i) = mod(i,2)*2 - 1;
%     
% %   Locations randomly in z = 0.5 plane
%     l = rand;
%     r = rand;
%     x0(:,i) = [l,r,.5]';
% 
% end

% % Ordered grid of xs in x-y plane
% l1 = linspace(0,1,pts1+1);
% it = 0;
% for i = 1:pts1
%     for j = 1:pts1
%         it = it+1;
%         x0(:,it) = [l1(i),l1(j),0.5]';
%     end
% end

% % Sphere xs
% p=4;
% nt = (p+1);
% [xs,wg] = lgwt(nt,-1,1);
% tht = acos(xs);
% np = 2*(p+1);
% pts = nt*np;
% dphi = 2*pi/np;
% phi = 0:dphi:dphi*(np-1)';
% it = 0;
% x0 = zeros(3,nt*np);
% f = x0;
% for i = 1:nt
%     for j = 1:np
%         it = it+1;
%         f(1,it) = (mod(it,2)*2 - 1)*sin(tht(i))*dphi*wg(i)*(0.25^2);
%         x0(1,it) = sin(tht(i))*cos(phi(j))/4+.5;
%         x0(2,it) = sin(tht(i))*sin(phi(j))/4+.5;
%         x0(3,it) = cos(tht(i))/4+.5;
%     end
% end
f = [   33.2417650606089       2.166587417781470E-008  -37.4980886580110     
   16.6208825034355        28.7882130043248       -37.4980886505129
  -16.6208825168548        28.7882129847437       -37.4980886479890
  -33.2417650399648       1.947638690896394E-008  -37.4980886527526
  -16.6208825446454       -28.7882129964959       -37.4980886601497
   16.6208825356381       -28.7882130138875       -37.4980886628842     
   84.0957521192305      -8.838722499782440E-009 -2.278202582629557E-007
   42.0478759789737        72.8290575286012      -2.244244899196829E-007
  -42.0478759286095        72.8290574519135      -2.270889122247098E-007
  -84.0957519359090      -1.432285365036645E-008 -2.323068415680900E-007
  -42.0478760430683       -72.8290576850160      -2.352984469127021E-007
   42.0478761044880       -72.8290577562195      -2.334762859125947E-007
   33.2417644466260       1.842359526230166E-008   37.4980889105340
   16.6208821957800        28.7882124705650        37.4980889056618
  -16.6208822118971        28.7882124510056        37.4980889011538     
  -33.2417644287216       1.627749617542111E-008   37.4980889017287
  -16.6208822397297       -28.7882124691772        37.4980889067020
   16.6208822279407       -28.7882124865906        37.4980889109994]';
x0=[  0.855153715697406       2.820947931273340E-003  0.401592096385261
  0.428987331773979       0.740962777300408       0.401592096454118
 -0.423345436021952       0.740962777271875       0.401592096463062     
 -0.849511819888657       2.820947875422286E-003  0.401592096402702
 -0.423345435962325      -0.735320881487066       0.401592096333578
  0.428987331827806      -0.735320881459748       0.401592096325081
   1.35047738261687       2.820947945247638E-003  2.820947119646463E-003
  0.676649165229340        1.16992565593949       2.820947157229190E-003
 -0.671007269472530        1.16992565589353       2.820947133082472E-003
  -1.34483548677237       2.820947856351632E-003  2.820947070234439E-003
 -0.671007269377574       -1.16428376012128       2.820947031983802E-003
  0.676649165309798       -1.16428376007835       2.820947057249106E-003
  0.855153714937626       2.820947975452442E-003 -0.395950199201968     
  0.428987331405093       0.740962776669471      -0.395950199223131
 -0.423345435611561       0.740962776640567      -0.395950199262878
 -0.849511819089884       2.820947918857954E-003 -0.395950199281911
 -0.423345435554447      -0.735320880768515      -0.395950199261015
  0.428987331456408      -0.735320880740825      -0.395950199220820  ]';
pts = numel(x0)/3;
% x=x0(:,1);
% Two points
pts = 2;
x0 = [.25,.5,.5
      .75,.5,.5]';
f= [1,0,0
   -1,0,0]';
% % One point
% pts = 1;
% f = zeros(3,pts);
% x0 = zeros(3,pts);
% x0(:,1) = [.75,.5,.5];
% f(:,1)= [1,0,0];


G = zeros(3,3,pts);
eps1 = 10;%1.8;%pi^0.5/tau^(1/3);

% Only calculating the real space part, as the Fourier part is done
% via fast ewald summation
for ps = 1:pts
for i = -bxs:bxs
    for j = -bxs:bxs
        for k = -bxs:bxs
            
            rc = i*bv(:,1) + j*bv(:,2) + k*bv(:,3);
            xn = (x - x0(:,ps)) + rc;
            
            G(:,:,ps) = G(:,:,ps) + HshortG(xn, eps1);
        end
    end
end
end

% ur = G(:,:,1)*f(:,1) + G(:,:,2)*f(:,2);

%% Doing the Fourier part

% First, make grid over periodic cell
np = 2^8;
h = 1/np;
% Basis vectors for our uniform grid
bv_gr = bv/(np-1);
bv_grx = bv_gr(:,1);
bv_gry = bv_gr(:,2);
bv_grz = bv_gr(:,3);

% Note: these are the basis vector differences in fractional space, which
% we integrate over since it's orthonormal and over cube [1,1,1];
hx = bv(1,1)/(np-1);
hy = hx;
hz = hy;

xs = zeros(np,np);
ys = zeros(np,np);

% Parameters dealing with shape of Gaussian interpolator
P = 80;%np;
m = 0.9*sqrt(pi*P);
w = P*hx/1.5;
cutpar = P*hx;
cut = 1e-16;
eta = -cutpar^2/log(cut)*eps1^2*2;%eps1^2*.35;%.35;%8.189752984211673E-002;%
par = (2*eps1^2/pi/eta)^1.5;
parexp = -2*eps1^2/eta;

% grid (for plots)z
for i = 1:np
    for j = 1:np
        rloc = bv_grx*(i-1) + bv_gry*(j-1);
        xs(i,j) = rloc(1);
        ys(i,j) = rloc(2);
    end
end        

H = zeros(3,np,np,np);
Hh = H;
Hht = H;
Ht = H;
Htint = H;
Hg = zeros(np,np,np);

supp_mat = zeros(3,1);
% Construct matrix containing all the combinations of points that fall
% within support. We want the term in the reduction from the exponential to
% be equal to cut
cut = 1e-16;
cutpar = sqrt(-log(cut)*eta/eps1^2/2);

% Due to the distortion of the grid, the minimum image convention can
% sometimes fail, choosing a point when a certain periodic image is closer.
% To fix this, we need to make sure that the cutoff distance is less than
% the minimum distance to the nearest image where this can fail, so that if
% it were to fail, the point is far enough away that it would've been
% truncated before then anyway. This distance is the shortest 1/2 distance
% between opposite faces, of which there are 3 pairs (easy to see in 2D,
% because the shortet distance parallel basis vectors (i.e. normal to both
% faces)

% Note that this really shouldn't ever be a problem because of the
% exponential decay. We should easily be able to have a cutoff below

% To find these distances, find the unit normal for a face using the two 
% basis vectors that dfine that face, and the distance will be the third
% vector dotted with the unit normal
n = cross(bvx,bvy);
r = abs(dot(n/norm(n),bvz));
n = cross(bvx,bvz);
r = min([r,abs(dot(n/norm(n),bvy))]);
n = cross(bvz,bvy);
r = min([r,abs(dot(n/norm(n),bvx))]);
r = r/2;

cutpar = 0.75;%min([r,cutpar]);

it = 0;
% probably slightly overkill, don't think I need to go over EVERY point
% Find indices in +- direction that are within cutoff (used for loops)
% Can do for the k too I suppose, but i only do the k loop once so it
% shouldn't actually change anything

% O(np^3)
for i = -floor(np/2):floor(np/2)
    for j = -floor(np/2):floor(np/2)
        for k = -floor(np/2):floor(np/2)
            rr = bv_grx*(i) + bv_gry*(j) + bv_grz*(k);
            r = norm(rr);
            if r<cutpar
                it = it + 1;
                supp_mat(:,it) = [i;j;k];
            end
            
        end
    end
end
supp_tot = it;

%% Actually starting to do the evaluating

% O(N*P^3)
for ip = 1:pts
% Find the grid point, rounded down
curijk = floor(bv_gr\x0(:,ip));

for it = 1:supp_tot
    indices = curijk + supp_mat(:,it)-1;
%   Managing periodicity
    if indices(1)<1
        iper = indices(1)+np;
    elseif indices(1)>np
        iper = indices(1)-np;
    else
        iper = indices(1);
    end
    
    if indices(2)<1
        jper = indices(2)+np;
    elseif indices(2)>np
        jper = indices(2)-np;
    else
        jper = indices(2);
    end

    if indices(3)<1
        kper = indices(3)+np;
    elseif indices(3)>np
        kper = indices(3)-np;
    else
        kper = indices(3);
    end
            
    indices=indices-1;
%   Use indices instead of iper etc, bc it will give distance to nearest
%   periodic image
    rr = x0(:,ip) ...
       - (bv_grx*indices(1) + bv_gry*indices(2) + bv_grz*indices(3));
    
    r2 = rr(1)^2 + rr(2)^2 + rr(3)^2;

    H(:,iper,jper,kper) = H(:,iper,jper,kper) + f(:,ip)*par*exp(parexp*r2);
    Hg(iper,jper,kper) = Hg(iper,jper,kper) + par*exp(parexp*r2);
%   Above is the same for the double layer
    if indices(1)==0 && indices(2)==5 && indices(3)==5
        l=1;
    end
end
end
surf(xs,ys,(real(squeeze(H(1,:,:,1+floor(0.5/bv_gr(3,3)))))))

Hh(1,:,:,:) = fftshift(fftn(H(1,:,:,:)));
Hh(2,:,:,:) = fftshift(fftn(H(2,:,:,:)));
Hh(3,:,:,:) = fftshift(fftn(H(3,:,:,:)));

HhtT = zeros(3,3,np,np,np);
HtT = HhtT;
HtintT = HhtT;

Hgh = fftshift(fftn(Hg));
Hght = zeros(3,3,np,np,np);
Hgt = zeros(3,3,np,np,np);
Hgtint=  Hgt;
aaa=0;
% Assuming even
% O(np^3)
for iw = -(np/2):np/2
    for jw = -(np/2):np/2
        for kw = -(np/2):np/2
            i = iw + np/2+1;
            j = jw + np/2+1;
            k = kw + np/2+1;
            
            k3 = kvx*iw + kvy*jw + kvz*kw;
            kn = norm(k3);
            if kn==0
                continue
            end
            if kn>40
                continue
            end
            
            B = 8*pi*(1+kn^2/4/eps1^2)/kn^4*(kn^2*eye(3) - k3*k3');
            B = B*exp((eta-1)*kn^2/4/eps1^2);
            aaa= max([aaa,max(max(abs(B)))]);
            Hht(1,i,j,k) = B(1,1)*Hh(1,i,j,k) + B(1,2)*Hh(2,i,j,k) + B(1,3)*Hh(3,i,j,k);
            Hht(2,i,j,k) = B(2,1)*Hh(1,i,j,k) + B(2,2)*Hh(2,i,j,k) + B(2,3)*Hh(3,i,j,k);
            Hht(3,i,j,k) = B(3,1)*Hh(1,i,j,k) + B(3,2)*Hh(2,i,j,k) + B(3,3)*Hh(3,i,j,k);
            
%           Pure Gij Green's function
            Hght(:,:,i,j,k) = B*Hgh(i,j,k);
%           Double layer amplification
            B = BT(k3, eps1)*exp((eta-1)*kn^2/4/eps1^2);
            HhtT(:,:,i,j,k) = B(:,:,1)*Hh(1,i,j,k) ...
                            + B(:,:,2)*Hh(2,i,j,k) ...
                            + B(:,:,3)*Hh(3,i,j,k);
        end
    end
end

% Single
Ht(1,:,:,:) = ifftn(ifftshift(Hht(1,:,:,:)));
Ht(2,:,:,:) = ifftn(ifftshift(Hht(2,:,:,:)));
Ht(3,:,:,:) = ifftn(ifftshift(Hht(3,:,:,:)));

% Double & single Greens
for  i = 1:3
    for j = 1:3
        HtT(i,j,:,:,:) = ifftn(ifftshift(HhtT(i,j,:,:,:)));
        Hgt(i,j,:,:,:) = ifftn(ifftshift(Hght(i,j,:,:,:)));
    end
end

% Find points with support around eval pt and only add in those
curijk = floor(bv_gr\x);

% O(N_x*P^3)
uu = zeros(3,1);
for it = 1:supp_tot
    indices = curijk + supp_mat(:,it);
%   Managing periodicity
    if indices(1)<1
        iper = indices(1)+np;
    elseif indices(1)>np
        iper = indices(1)-np;
    else
        iper = indices(1);
    end
    
    if indices(2)<1
        jper = indices(2)+np;
    elseif indices(2)>np
        jper = indices(2)-np;
    else
        jper = indices(2);
    end

    if indices(3)<1
        kper = indices(3)+np;
    elseif indices(3)>np
        kper = indices(3)-np;
    else
        kper = indices(3);
    end
    indices=indices-1;
            
    rr = x ...
       - (bv_grx*indices(1) + bv_gry*indices(2) + bv_grz*indices(3));
   
    r2 = rr(1)^2 + rr(2)^2 + rr(3)^2;
    Htint(:,iper,jper,kper) = Ht(:,iper,jper,kper)*par*exp(parexp*r2);
    
%   Double
    HtintT(:,:,iper,jper,kper) = HtT(:,:,iper,jper,kper)*par*exp(parexp*r2);
    
%   Single
    Hgtint(:,:,iper,jper,kper) = Hgt(:,:,iper,jper,kper)*par*exp(parexp*r2);
end

% Jacobian of transformation from fractional to real space (transfo is just
% the basis vectors x = bv*u).
J = det(bv);

uu = zeros(3,1);
uu(1) = sum(sum(sum(Htint(1,:,:,:))))*h*h*h*J;
uu(2) = sum(sum(sum(Htint(2,:,:,:))))*h*h*h*J;
uu(3) = sum(sum(sum(Htint(3,:,:,:))))*h*h*h*J;
disp(real(uu))

Dij = zeros(3,3);

for i = 1:3
    for j = 1:3
        Dij(i,j) = sum(sum(sum(HtintT(i,j,:,:,:))))*h*h*h*J;
    end
end
% disp(Dij)

Sij = zeros(3,3);
for i = 1:3
    for j = 1:3
        Sij(i,j) = sum(sum(sum(Hgtint(i,j,:,:,:))))*h*h*h*J;
    end
end
% disp(Sij)



%% old stuff

% We need to find the nearest periodic image of each grid point to each
% point force. Just construct the grid ahead of time here
% gridnn = zeros(3,pts,np,np,np);
% evalnn = zeros(3,1,np,np,np);

% Note that for very distorted grids, this could cause us to exclude points
% that are actually nearer... See solution below in supp_mat

% % Loop over grid points
% for i = 0:np-1
%     for j = 0:np-1
%         for k = 0:np-1
%             tmpr = bv_grx*i + bv_gry*j + bv_grz*k;
% %           Loop over point forces
%             for  ip = 1:pts
%                 rr = x0(:,ip) - tmpr;
% %               If the distance in 1D is greater than half a basis vector,
% %               the closest image will be in a different box.
% %               Start with cubic for now
%                 gridnn(:,ip,i+1,j+1,k+1) = round(-rr); % only for cubic, 1x1x1
%             end
%         end
%     end
% end

% % Also need to do this for the evaluation points, although these should
% % end up being the same as force locations in actual implementation
% for i = 0:np-1
%     for j = 0:np-1
%         for k = 0:np-1
%             tmpr = bv_grx*i + bv_gry*j + bv_grz*k;
%             rr = x - tmpr;
% %           If the distance in 1D is greater than half a basis vector,
% %           the closest image will be in a different box.
% %           Start with cubic for now
%             evalnn(:,ip,i+1,j+1,k+1) = round(-rr); % only for cubic, 1x1x1
%         end
%     end
% end


%% Functions
% Double layer amplification factor
function B = BT(k, eps1)
kn = norm(k);
D = zeros(3,3,3);
D(1,1,1) = -2/(kn^4)*k(1)*k(1)*k(1) + 1/(kn^2)*(k(1) + k(1) + k(1));
D(2,1,1) = -2/(kn^4)*k(2)*k(1)*k(1) + 1/(kn^2)*(       k(2)       );
D(3,1,1) = -2/(kn^4)*k(3)*k(1)*k(1) + 1/(kn^2)*(       k(3)       );
D(1,2,1) = -2/(kn^4)*k(1)*k(2)*k(1) + 1/(kn^2)*(              k(2));
D(2,2,1) = -2/(kn^4)*k(2)*k(2)*k(1) + 1/(kn^2)*(k(1)              );
D(3,2,1) = -2/(kn^4)*k(3)*k(2)*k(1);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(1,3,1) = -2/(kn^4)*k(1)*k(3)*k(1) + 1/(kn^2)*(              k(3));
D(2,3,1) = -2/(kn^4)*k(2)*k(3)*k(1);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(3,3,1) = -2/(kn^4)*k(3)*k(3)*k(1) + 1/(kn^2)*(k(1)              );
D(1,1,2) = -2/(kn^4)*k(1)*k(1)*k(2) + 1/(kn^2)*(k(2)              );
D(2,1,2) = -2/(kn^4)*k(2)*k(1)*k(2) + 1/(kn^2)*(              k(1));
D(3,1,2) = -2/(kn^4)*k(3)*k(1)*k(2);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(1,2,2) = -2/(kn^4)*k(1)*k(2)*k(2) + 1/(kn^2)*(       k(1)       );
D(2,2,2) = -2/(kn^4)*k(2)*k(2)*k(2) + 1/(kn^2)*(k(2) + k(2) + k(2));
D(3,2,2) = -2/(kn^4)*k(3)*k(2)*k(2) + 1/(kn^2)*(       k(3)       );
D(1,3,2) = -2/(kn^4)*k(1)*k(3)*k(2);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(2,3,2) = -2/(kn^4)*k(2)*k(3)*k(2) + 1/(kn^2)*(              k(3));
D(3,3,2) = -2/(kn^4)*k(3)*k(3)*k(2) + 1/(kn^2)*(k(2)              );
D(1,1,3) = -2/(kn^4)*k(1)*k(1)*k(3) + 1/(kn^2)*(k(3)              );
D(2,1,3) = -2/(kn^4)*k(2)*k(1)*k(3);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(3,1,3) = -2/(kn^4)*k(3)*k(1)*k(3) + 1/(kn^2)*(              k(1));
D(1,2,3) = -2/(kn^4)*k(1)*k(2)*k(3);%+1/(kn^2)*(k(1) + k(1) + k(1));
D(2,2,3) = -2/(kn^4)*k(2)*k(2)*k(3) + 1/(kn^2)*(k(3)              );
D(3,2,3) = -2/(kn^4)*k(3)*k(2)*k(3) + 1/(kn^2)*(              k(2));
D(1,3,3) = -2/(kn^4)*k(1)*k(3)*k(3) + 1/(kn^2)*(       k(1)       );
D(2,3,3) = -2/(kn^4)*k(2)*k(3)*k(3) + 1/(kn^2)*(       k(2)       );
D(3,3,3) = -2/(kn^4)*k(3)*k(3)*k(3) + 1/(kn^2)*(k(3) + k(3) + k(3));
B1 = 8*pi*(1 + kn^2/eps1^2*0.25);
B = B1*D*1i;
end


% Hasimotos
function G = HshortG(x, eps1)
r = norm(x);
er = r*eps1;
C = erfc(er) - 2/sqrt(pi)*er*exp(-er^2);
D = erfc(er) + 2/sqrt(pi)*er*exp(-er^2);
G = eye(3)*C/r + x*x'*D/r^3;
end

% function G = HspecG(k, eps1)
% kn = norm(k);
% w = kn/eps1;
% G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2      )*(kn^2*eye(3) - k*k')*exp(-w^2/4);
% end
% 
% function B = partHspec(k,eps1)
% kn = norm(k);
% B = 8*pi*(1 + (kn/eps1)^2)/(kn^4)*(kn^2*eye(3) - k*k');
% end

%%% See Greengard and Lee (via tornberg) for speeding up

% First, define the lattice (columns vectors)
bv = [1 0 0; 0 1 0; 0 0 1];
bvx = bv(:,1);
bvy = bv(:,2);
bvz = bv(:,3);
x0 = [.25,.5,.5
      .75,.5,.5]';
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

eps1 = pi^0.5/tau^(1/3);

f= [1,0,0
   -1,0,0]';

pts = 2;
bxs = 2;
G = zeros(3,3,pts);

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

ur = G(:,:,1)*f(:,1) + G(:,:,2)*f(:,2);

%% Doing the Fourier part

% First, make grid over periodic cell
np = 2^6;
h = 1/np;
% Basis vectors for our uniform grid
bv_gr = bv*h;
bv_grx = bv_gr(:,1);
bv_gry = bv_gr(:,2);
bv_grz = bv_gr(:,3);

% Note: these are the basis vector differences in fractional space, which
% we integrate over since it's orthonormal and over cube [1,1,1];
hx = 1/np;
hy = hx;
hz = hy;

xs = zeros(np,np);
ys = zeros(np,np);

% Parameters dealing with shape of Gaussian interpolator
P = np;
m = 0.9*sqrt(pi*P);
w = P*h/1.5;
eta = (2*w*eps1/m)^2;
par = (2*eps1^2/pi/eta)^1.5;
parexp = -2*eps1^2/eta;

% grid (for plots)
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

supp_mat = zeros(3,1);
% Construct matrix containing all the combinations of points that fall
% within support. We want the term in the reduction from the exponential to
% be equal to cut
cut = 1e-12;
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

cutpar = min([r,cutpar]);

it = 0;
% probably slightly overkill, don't think I need to go over EVERY point
% Find indices in +- direction that are within cutoff (used for loops)
% Can do for the k too I suppose, but i only do the k loop once so it
% shouldn't actually change anything

% O(np^3)
for i = -floor(np/2):floor(np/2)
    for j = -floor(np/2):floor(np/2)
        for k = -floor(np/2):floor(np/2)
            rr = bv_grx*i + bv_gry*j + bv_grz*k;
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
            
%   Use indices instead of iper etc, bc it will give distance to nearest
%   periodic image
    rr = x0(:,ip) ...
       - (bv_grx*indices(1) + bv_gry*indices(2) + bv_grz*indices(3));
    
    r2 = rr(1)^2 + rr(2)^2 + rr(3)^2;

    H(:,iper,jper,kper) = H(:,iper,jper,kper) + f(:,ip)*par*exp(parexp*r2);
end
end

Hh(1,:,:,:) = fftshift(fftn(H(1,:,:,:)));
Hh(2,:,:,:) = fftshift(fftn(H(2,:,:,:)));
Hh(3,:,:,:) = fftshift(fftn(H(3,:,:,:)));

% Assuming even
% O(np^3)
for iw = -(np/2):np/2-1
    for jw = -(np/2):np/2-1
        for kw = -(np/2):np/2-1
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
            Hht(1,i,j,k) = B(1,1)*Hh(1,i,j,k) + B(1,2)*Hh(2,i,j,k) + B(1,3)*Hh(3,i,j,k);
            Hht(2,i,j,k) = B(2,1)*Hh(1,i,j,k) + B(2,2)*Hh(2,i,j,k) + B(2,3)*Hh(3,i,j,k);
            Hht(3,i,j,k) = B(3,1)*Hh(1,i,j,k) + B(3,2)*Hh(2,i,j,k) + B(3,3)*Hh(3,i,j,k);
        end
    end
end

Ht(1,:,:,:) = ifftn(ifftshift(Hht(1,:,:,:)));
Ht(2,:,:,:) = ifftn(ifftshift(Hht(2,:,:,:)));
Ht(3,:,:,:) = ifftn(ifftshift(Hht(3,:,:,:)));


% Find points with support around eval pt and only add in those
curijk = floor(bv_gr\x);

% O(N_x*P^3)
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
            
    rr = x ...
       - (bv_grx*indices(1) + bv_gry*indices(2) + bv_grz*indices(3));
   
    r2 = rr(1)^2 + rr(2)^2 + rr(3)^2;
    Htint(:,iper,jper,kper) = Ht(:,iper,jper,kper)*par*exp(parexp*r2);
end

% Jacobian of transformation from fractional to real space (transfo is just
% the basis vectors x = bv*u).
J = det(bv);

uu = zeros(3,1);
uu(1) = sum(sum(sum(Htint(1,:,:,:))))*hx*hy*hz*J;
uu(2) = sum(sum(sum(Htint(2,:,:,:))))*hx*hy*hz*J;
uu(3) = sum(sum(sum(Htint(3,:,:,:))))*hx*hy*hz*J;
disp(uu)

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
% Hasimotos
function G = HshortG(x, eps1)
r = norm(x);
er = r*eps1;
C = erfc(er) - 2/sqrt(pi)*er*exp(-er^2);
D = erfc(er) + 2/sqrt(pi)*er*exp(-er^2);
G = eye(3)*C/r + x*x'*D/r^3;
end

function G = HspecG(k, eps1)
kn = norm(k);
w = kn/eps1;
G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2      )*(kn^2*eye(3) - k*k')*exp(-w^2/4);
end

function B = partHspec(k,eps1)
kn = norm(k);
B = 8*pi*(1 + (kn/eps1)^2)/(kn^4)*(kn^2*eye(3) - k*k');
end

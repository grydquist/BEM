% Surface constants
xmn1 = x1c;
xmn2 = x2c;
xmn3 = x3c;
p = size(xmn1);
p = sqrt(p(1)) - 1;

% Constants for phi
np = 2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';

% Constants for theta
nt = (p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);
% Harmonics evaluated at the given thetas/phis
Yt = SpHarmTNew(p,th,ph);

% Making the fine grid
q= 2*p;
ntf = (q+1);
npf = 2*(q+1);

dphif = 2*pi/npf;
phif = 0:dphif:dphif*(npf-1)';
[xsf,wgf] = lgwt(ntf,-1,1);
thtf = acos(xsf);
[phf,thf] = meshgrid(phif,thtf);
Ytf = SpHarmTNew(p,thf,phf);

% Get all the derivatives of the Spherical Harmonics up front
Ytd1 = zeros(ntf,npf,(p+1)^2);

it = 0;
im = 0;
for n = 0:p
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
%       Theta derivative
        Ytd1(:,:,it) = ThetDer(Ytf,phf,n,m,1);
    end
end

% Construct Grids
x1 = real(SpHReconst(x1c,Yt));
x2 = real(SpHReconst(x2c,Yt));
x3 = real(SpHReconst(x3c,Yt));
x1f = real(SpHReconst(x1c,Ytf));
x2f = real(SpHReconst(x2c,Ytf));
x3f = real(SpHReconst(x3c,Ytf));

% Total derivatives, in Cartesian
dxp = zeros(3,ntf,npf);
dxt = dxp;

% Calcaulating surface derivatives
% Loop over harmonics
ih=0;
it=0;

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
        f1 = xmn1(it);
        f2 = xmn2(it);
        f3 = xmn3(it);
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


%% Construct attenuation coefficients (perfect low pass filter)
amn = xmn1;
amn(:) = 1;

% Filter at n = p/3
ncut = ceil(p/4);
amn(1:(ncut+1)^2) = 0;

% Calculate quality measure
it = 1;
E = 0;
for n=1:p
    for m = -n:n
        xcur = [xmn1(it), xmn2(it), xmn3(it)];
        E = E + amn(it)*norm(xcur)^2;
        it = it + 1;
    end
end

% First, let's evaluate the adaptive cutoff, according to Sorgentone:

N2 = zeros(1,p);
for no = 1:p
    N2(no) = 0;
    it = no^2;
    for n = no:p
        for m = -n:n
            it = it+1;
            N2(no) = N2(no) + norm([xmn1(it);xmn2(it);xmn3(it)]);
        end
    end
end
N2p = N2/N2(1);
ncut = N2p( 1:find( N2p < 0.2, 1 ) );
ncut = numel(ncut) - 1;
amn(:) = 1;
amn(1:(ncut)^2) = 0;

%% Now loop through all the points and advance them according to the optimality measure

dtau = .25;

for loop=1:5
for i = 1:ntf
    for j = 1:npf
%       Need the normal vector
%%%%%%%%%% This normal vector needs updating!!
        nk = cross(dxt(:,i,j),dxp(:,i,j)); 
        nk = -nk./norm(nk);
        nk = real(nk);
       
       
%       Gradient of energy
        delE = [0;0;0];
        for n = ncut:p
            ih = n + 1;
            it = (n)^2;
            im=0;
            Ypcur = Ytf{ih};
            for m = -n:n
                im = im+1;
                it = it + 1;
                Ym = squeeze(Ypcur(im,:,:));
                delE(1) = delE(1) + xmn1(it)*Ym(i,j);
                delE(2) = delE(2) + xmn2(it)*Ym(i,j);
                delE(3) = delE(3) + xmn3(it)*Ym(i,j);
            end
        end

        ycur = [x1f(i,j);x2f(i,j);x3f(i,j)] - dtau*(eye(3) - nk*nk')*real(delE);
        x1f(i,j) = ycur(1);
        x2f(i,j) = ycur(2);
        x3f(i,j) = ycur(3);
       
    end
end
xmn1 = SpT(Ytf,x1f,thtf,phif,p);
xmn2 = SpT(Ytf,x2f,thtf,phif,p);
xmn3 = SpT(Ytf,x3f,thtf,phif,p);

end




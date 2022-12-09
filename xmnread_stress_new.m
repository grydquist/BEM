%% Pre allocation
fclose all;
% Harmonics evaluated at fine grid
ntf = 25;
npf = 50;
thtf = linspace(.001,pi-.001,ntf);
% [xsf,wgf] = lgwt(ntf,-1,1);
% thtf = acos(xsf);
phif = linspace(0,2*pi,npf);
[phf,thf] = meshgrid(phif,thtf);

dxp = zeros(3,ntf,npf);
xf = dxp;
dxt = dxp;
dxp2= dxp;
dxt2= dxp;
dxtp= dxp;
dxp3 = dxp;
dxt3 = dxp;
dxt2p= dxp;
dxtp2= dxp;
dxp4 = dxp;
dxtp3 = dxp;
dxt2p2 = dxp;
dxt3p = dxp;
dxt4 = dxp;

% Normal Vector (inward)
nk = zeros(3,ntf,npf);

J = zeros(ntf,npf);
JR = J;
ps1 = J;
ps2 = J;
Ar = J;
Sh = J;
taup1 = J;
taup2 = J;
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
myf = zeros(3,ntf,npf);
fab = myf;

% Read Param file to get simulation info
dir = 'fortran/dat/p12S250H45/cell_2/';
% dir = 'fortran/dat/TEST/cell_1/';
% dir = 'fortran/dat/IndNewSR250/cell_1/';
% dir = 'fortran/dat/TurbRes/Turb20/cell_1/';
fid = fopen(strcat(dir,'Params'));
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    switch tline
        case 'p'
            p = str2double(fgetl(fid));
            p = floor(p);
        case 'dt'
            ts = str2double(fgetl(fid));
        case 'dt_inc'
            dt_inc = str2double(fgetl(fid));
            dt_inc = floor(dt_inc);
        case 'Ed'
            Ed = str2double(fgetl(fid));
        case 'Ca'
            Ca = str2double(fgetl(fid));
        case 'c0'
            c0 = str2double(fgetl(fid));
        case 'lambda'
            lam = str2double(fgetl(fid));
        case 'Eb'
            Ebs = str2double(fgetl(fid));
    end
    
    tline = fgetl(fid);
end
fclose(fid);

B = 2/Ca;
C = Ed*B;
Eb = 2*B*Ebs;

%% Actual Reading
% Reads the txt file output from the fortran code for shape and calculates
% stresses and whatnot

% Get total timesteps outputted
fid = fopen(strcat(dir,'maxdt'));
tts = str2double(fgetl(fid));
fclose(fid);
tts = floor(tts);

% How many timesteps to skip
tincr = .5;

incr = floor(tincr/ts);
% Round down to fit w/ dt_inc
incr = incr - mod(incr,dt_inc);
if(incr == 0); incr = dt_inc; end
tsteps = floor(tts/incr) + 1;

% Time
t = zeros(tsteps,1);

% Material point to track
xtop = zeros(tsteps,1);
ytop = zeros(tsteps,1);
ztop = zeros(tsteps,1);
Dij = ytop;
incl = Dij;
elxa = zeros(tsteps,202);
elza = elxa;

% Number total values in a time step
lent1 = 6*(p+1)^2;
% length of 1 collection in one time step
tot = (p+1)^2;

% Order of spherical harmonics
Ex = zeros(p+1,tsteps);
Eu = zeros(p+1,tsteps);

% Evaluation of spherical harmonics for interpolation
tmpt = linspace(0,pi,101);tmpp = linspace(0,2*pi,101);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmTNew(p,tmpth,tmpph);

tq = linspace(0,pi,15);pq = tq;
[ttq,ppq] = meshgrid(tq,pq);
Yqv = SpHarmTNew(p,ttq,ppq);
Ytf = SpHarmTNew(p,thf,phf);
% Spherical harmonic evaluated at top of sphere
Ytrc = SpHarmTNew(p,0,0);

% Max areal strain
Jmax = t;
Jmin = t;
pst1 = t;
pst2 = t;
Shmax= t;
Shmin= t;

% Get all the derivatives of the Spherical Harmonics up front
Ytd1 = zeros(ntf,npf,(p+1)^2);
Ytd2 = Ytd1;
Ytd3 = Ytd1;
Ytd4 = Ytd1;

it = 0;
im = 0;
for n = 0:p
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
        Ytd1(:,:,it) = ThetDer(Ytf,phf,n,m,1);
        Ytd2(:,:,it) = ThetDer(Ytf,phf,n,m,2);
        Ytd3(:,:,it) = ThetDer(Ytf,phf,n,m,3);
        Ytd4(:,:,it) = ThetDer(Ytf,phf,n,m,4);
    end
end

%% Actual plotting
disp('Start!')
% Do timesteps
for i = 1:incr:tts + 1
% Current time
    t((i-1)/incr + 1) = i*ts - ts;
    
%   Read in data file of current timestep, making exception for the first
    if(i == 1)
        file = 'x_00001';
    else
        file = 'x_';
        for fl = 1:(4-floor(log10(i-1)))
            file = strcat(file,'0');
        end
        file = strcat(file,num2str(i-1));
    end
    
%   All the data in the ts
    fID = fopen(strcat(dir,file));
    raw = fscanf(fID,'%f');
    fclose(fID);
    
%   Individual directional coefficients, not distinct between real/imag)
    x1t = raw(1:3:end);
    x2t = raw(2:3:end);
    x3t = raw(3:3:end);
    
%   Real and imaginary separation
    x1c = x1t(1:tot) + x1t(tot+1:end)*1i;
    x2c = x2t(1:tot) + x2t(tot+1:end)*1i;
    x3c = x3t(1:tot) + x3t(tot+1:end)*1i;
    
%%  Actual calculation of things like area element, stresses, etc after here
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

%%  Calcaulating surface derivatives
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
            f1 = x1c(it);
            f2 = x2c(it);
            f3 = x3c(it);
    %       Current values of harmonic at degree and order
            Ym = squeeze(Ypcur(im,:,:));
    %       Theta derivative
            td1 = Ytd1(:,:,it);

    %       Distance from zero
            xf(1,:,:) = squeeze(xf(1,:,:)) + f1*Ym;
            xf(2,:,:) = squeeze(xf(2,:,:)) + f2*Ym;
            xf(3,:,:) = squeeze(xf(3,:,:)) + f3*Ym;
            xf = real(xf);

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


    %       Fourth order derivatives (for shear force derivatives)
    %       Fourth phi derivative
            dxp4(1,:,:) = squeeze(dxp4(1,:,:)) + f1*m^4*Ym;
            dxp4(2,:,:) = squeeze(dxp4(2,:,:)) + f2*m^4*Ym;
            dxp4(3,:,:) = squeeze(dxp4(3,:,:)) + f3*m^4*Ym;

    %       Mixed derivatives
            dxtp3(1,:,:) = squeeze(dxtp3(1,:,:)) + f1*-1i*m^3*td1;
            dxtp3(2,:,:) = squeeze(dxtp3(2,:,:)) + f2*-1i*m^3*td1;
            dxtp3(3,:,:) = squeeze(dxtp3(3,:,:)) + f3*-1i*m^3*td1;

            dxt2p2(1,:,:) = squeeze(dxt2p2(1,:,:)) + f1*-m^2*td2;
            dxt2p2(2,:,:) = squeeze(dxt2p2(2,:,:)) + f2*-m^2*td2;
            dxt2p2(3,:,:) = squeeze(dxt2p2(3,:,:)) + f3*-m^2*td2;

            dxt3p(1,:,:) = squeeze(dxt3p(1,:,:)) + f1*1i*m*td3;
            dxt3p(2,:,:) = squeeze(dxt3p(2,:,:)) + f2*1i*m*td3;
            dxt3p(3,:,:) = squeeze(dxt3p(3,:,:)) + f3*1i*m*td3;

    %       Fourth theta derivative
            td4 = Ytd4(:,:,it);
            dxt4(1,:,:) = squeeze(dxt4(1,:,:)) + f1*td4;
            dxt4(2,:,:) = squeeze(dxt4(2,:,:)) + f2*td4;
            dxt4(3,:,:) = squeeze(dxt4(3,:,:)) + f3*td4;
        end
    end
    % These shouldn't be imaginary, but some creeps in        
    dxt   = real(dxt);
    dxt2  = real(dxt2);
    dxt3  = real(dxt3);

    dxp   = real(dxp);
    dxp2  = real(dxp2);
    dxp3  = real(dxp3);

    dxtp  = real(dxtp);
    dxt2p = real(dxt2p);
    dxtp2 = real(dxtp2);

    dxp4 = real(dxp4);
    dxtp3 = real(dxtp3);
    dxt2p2 = real(dxt2p2);
    dxt3p = real(dxt3p);
    dxt4 = real(dxt4);
    
%% And now stresses

for i2 = 1:ntf
    for j = 1:npf
        
%       Normal vector (inward)
        nk(:,i2,j) = cross(dxt(:,i2,j),dxp(:,i2,j));
        nk(:,i2,j) = -nk(:,i2,j)./norm(nk(:,i2,j));
        nk(:,i2,j) = real(nk(:,i2,j));
        
%       Jacobian via fundamental forms
        E = dot(dxt(:,i2,j),dxt(:,i2,j));
        F = dot(dxt(:,i2,j),dxp(:,i2,j));
        G = dot(dxp(:,i2,j),dxp(:,i2,j));
        J(i2,j) = sqrt(E.*G - F.^2);

%       Calculate mean curvature with some fundamental forms (appendices in
%       Veeranpani paper)
        L = dot(dxt2(:,i2,j),-nk(:,i2,j));
        M = dot(dxtp(:,i2,j),-nk(:,i2,j));
        N = dot(dxp2(:,i2,j),-nk(:,i2,j));
%       Curvature numerator
        D = E*N - 2*F*M + G*L;
        k = 0.5*D/J(i2,j)^2;
%       Gaussian curvature
        kG = (L*N - M*M)/(J(i2,j)*J(i2,j));

%       First let's do all the cross products we need
        ct2_p = cross(dxt2(:,i2,j),dxp(:,i2,j));
        ct_tp = cross(dxt (:,i2,j),dxtp(:,i2,j));
        ctp_p = cross(dxtp(:,i2,j),dxp(:,i2,j));
        ct_p2 = cross(dxt (:,i2,j),dxp2(:,i2,j));
        
        ct3_p = cross(dxt3(:,i2,j),dxp(:,i2,j));
        ct2_tp = cross(dxt2(:,i2,j),dxtp(:,i2,j));
        ct_t2p = cross(dxt(:,i2,j),dxt2p(:,i2,j));
        
        ct2p_p = cross(dxt2p(:,i2,j),dxp(:,i2,j));
        ct2_p2 = cross(dxt2(:,i2,j),dxp2(:,i2,j));
        ct_tp2 = cross(dxt(:,i2,j),dxtp2(:,i2,j));
        
        ctp2_p = cross(dxtp2(:,i2,j),dxp(:,i2,j));
        ctp_p2 = cross(dxtp(:,i2,j),dxp2(:,i2,j));
        ct_p3 = cross( dxt(:,i2,j),dxp3(:,i2,j));
        
%       Jacobian Partials        
        Jt = dot(-nk(:,i2,j),ct2_p + ct_tp);
        Jp = dot(-nk(:,i2,j),ctp_p + ct_p2);

%       Normal vector partials outward
        dnt = (1/J(i2,j))*(ct2_p + ct_tp - Jt*-nk(:,i2,j));
        dnp = (1/J(i2,j))*(ctp_p + ct_p2 - Jp*-nk(:,i2,j));
        
%       Second Jacobian derivs
        Jt2 = dot(dnt,ct2_p + ct_tp) ...
            + dot(-nk(:,i2,j), ct3_p + 2*ct2_tp + ct_t2p);
        Jtp = dot(dnp,ct2_p + ct_tp) ...
            + dot(-nk(:,i2,j), ct2p_p +  ct2_p2 + ct_tp2);
        Jp2 = dot(dnp,ctp_p + ct_p2) ...
            + dot(-nk(:,i2,j), ctp2_p + 2*ctp_p2 + ct_p3);
        
%       Second normal derivatives        
        dnt2 = (1/J(i2,j))*(-2*Jt*dnt - -nk(:,i2,j)*Jt2 + ct3_p + 2*ct2_tp + ct_t2p);
        dnp2 = (1/J(i2,j))*(-2*Jp*dnp - -nk(:,i2,j)*Jp2 + ct_p3 + 2*ctp_p2 + ctp2_p);
        dntp = (1/J(i2,j))*(-Jp*dnt - Jt*dnp - -nk(:,i2,j)*Jtp + ct2p_p + ct2_p2 + ct_tp2);
     
%       Fundamental form partials        
        Et  = 2*dot(dxt2(:,i2,j),dxt(:,i2,j));
        Ep  = 2*dot(dxtp(:,i2,j),dxt(:,i2,j));
        Et2 = 2*(dot(dxt3(:,i2,j),dxt(:,i2,j)) + dot(dxt2(:,i2,j),dxt2(:,i2,j)));
        Ep2 = 2*(dot(dxtp2(:,i2,j),dxt(:,i2,j))+ dot(dxtp(:,i2,j),dxtp(:,i2,j)));
        Etp = 2*(dot(dxt2p(:,i2,j),dxt(:,i2,j))+ dot(dxt2(:,i2,j),dxtp(:,i2,j)));
        
        Gt  = 2*dot(dxtp(:,i2,j),dxp(:,i2,j));
        Gp  = 2*dot(dxp (:,i2,j),dxp2(:,i2,j));
        Gt2 = 2*(dot(dxt2p(:,i2,j),dxp(:,i2,j))+ dot(dxtp(:,i2,j),dxtp(:,i2,j)));
        Gp2 = 2*(dot(dxp3(:,i2,j),dxp(:,i2,j)) + dot(dxp2(:,i2,j),dxp2(:,i2,j)));
        Gtp = 2*(dot(dxtp2(:,i2,j),dxp(:,i2,j))+ dot(dxtp(:,i2,j),dxp2(:,i2,j)));
        
        Ft  = dot(dxt2(:,i2,j),dxp(:,i2,j)) + dot(dxt(:,i2,j),dxtp(:,i2,j));
        Fp  = dot(dxtp(:,i2,j),dxp(:,i2,j)) + dot(dxt(:,i2,j),dxp2(:,i2,j));
        Ft2 = dot(dxt3(:,i2,j),dxp(:,i2,j)) + 2*dot(dxt2(:,i2,j),dxtp(:,i2,j)) + dot(dxt(:,i2,j),dxt2p(:,i2,j));
        Fp2 = dot(dxtp2(:,i2,j),dxp(:,i2,j))+ 2*dot(dxtp(:,i2,j),dxp2(:,i2,j)) + dot(dxt(:,i2,j),dxp3(:,i2,j));
        Ftp = dot(dxt2p(:,i2,j),dxp(:,i2,j))+ dot(dxtp(:,i2,j),dxtp(:,i2,j)) + dot(dxt2(:,i2,j),dxp2(:,i2,j)) + dot(dxt(:,i2,j),dxtp2(:,i2,j));
        
        Lt  = dot(dxt3(:,i2,j),-nk(:,i2,j))  + dot(dxt2(:,i2,j),dnt);
        Lp  = dot(dxt2p(:,i2,j),-nk(:,i2,j)) + dot(dxt2(:,i2,j),dnp);
        Lt2 = dot(dxt4(:,i2,j),-nk(:,i2,j))  + 2*dot(dxt3(:,i2,j),dnt) + dot(dxt2(:,i2,j),dnt2);
        Lp2 = dot(dxt2p2(:,i2,j),-nk(:,i2,j))+ 2*dot(dxt2p(:,i2,j),dnp) + dot(dxt2(:,i2,j),dnp2);
        Ltp = dot(dxt3p(:,i2,j),-nk(:,i2,j)) + dot(dxt2p(:,i2,j),dnt) + dot(dxt3(:,i2,j),dnp) + dot(dxt2(:,i2,j),dntp);
        
        Nt  = dot(dxtp2(:,i2,j),-nk(:,i2,j)) + dot(dxp2(:,i2,j),dnt);
        Np  = dot(dxp3(:,i2,j),-nk(:,i2,j))  + dot(dxp2(:,i2,j),dnp);
        Nt2 = dot(dxt2p2(:,i2,j),-nk(:,i2,j))+ 2*dot(dxtp2(:,i2,j),dnt) + dot(dxp2(:,i2,j),dnt2);
        Np2 = dot(dxp4(:,i2,j),-nk(:,i2,j))  + 2*dot(dxp3(:,i2,j),dnp) + dot(dxp2(:,i2,j),dnp2);
        Ntp = dot(dxtp3(:,i2,j),-nk(:,i2,j)) + dot(dxtp2(:,i2,j),dnp) + dot(dxp3(:,i2,j),dnt) + dot(dxp2(:,i2,j),dntp);
        
        Mt  = dot(dxt2p(:,i2,j),-nk(:,i2,j)) + dot(dxtp(:,i2,j),dnt);
        Mp  = dot(dxtp2(:,i2,j),-nk(:,i2,j)) + dot(dxtp(:,i2,j),dnp);
        Mt2 = dot(dxt3p(:,i2,j),-nk(:,i2,j)) + 2*dot(dxt2p(:,i2,j),dnt) + dot(dxtp(:,i2,j),dnt2);
        Mp2 = dot(dxtp3(:,i2,j),-nk(:,i2,j)) + 2*dot(dxtp2(:,i2,j),dnp) + dot(dxtp(:,i2,j),dnp2);
        Mtp = dot(dxt2p2(:,i2,j),-nk(:,i2,j))+ dot(dxt2p(:,i2,j),dnp) + dot(dxtp2(:,i2,j),dnt) + dot(dxtp(:,i2,j),dntp);

%       Curvature derivatives
%       Numerator derivatives
        Dt = Et*N + E*Nt - 2*Ft*M - 2*F*Mt + Gt*L + G*Lt;
        Dp = Ep*N + E*Np - 2*Fp*M - 2*F*Mp + Gp*L + G*Lp;
        
        Dt2 = Et2*N + 2*Et*Nt + E*Nt2 - 2*Ft2*M - 4*Ft*Mt - 2*F*Mt2 + Gt2*L + 2*Gt*Lt + G*Lt2;
        Dp2 = Ep2*N + 2*Ep*Np + E*Np2 - 2*Fp2*M - 4*Fp*Mp - 2*F*Mp2 + Gp2*L + 2*Gp*Lp + G*Lp2;
        Dtp = Etp*N + Et*Np + Ep*Nt + E*Ntp - 2*Ftp*M - 2*Ft*Mp - 2*Fp*Mt - 2*F*Mtp + Gtp*L + Gt*Lp + Gp*Lt + G*Ltp;
        
        kt = 0.5*Dt*J(i2,j)^-2 - 2*k*Jt/J(i2,j);
        kp = 0.5*Dp*J(i2,j)^-2 - 2*k*Jp/J(i2,j);
            
        kt2 = 0.5*Dt2*J(i2,j)^-2 - 2*Jt*J(i2,j)^-3*Dt +3*J(i2,j)^-4*Jt^2*D - J(i2,j)^-3*Jt2*D;
        kp2 = 0.5*Dp2*J(i2,j)^-2 - 2*Jp*J(i2,j)^-3*Dp +3*J(i2,j)^-4*Jp^2*D - J(i2,j)^-3*Jp2*D;
        ktp = 0.5*Dtp*J(i2,j)^-2 - Jt*J(i2,j)^-3*Dp - Jp*J(i2,j)^-3*Dt +3*J(i2,j)^-4*Jt*Jp*D - J(i2,j)^-3*Jtp*D;
        
%       Metric tensor
        g = [E,F;F,G];
%       Contravariant g
        gn = [G,-F;-F,E]/(J(i2,j)^2);
%       Partials of COTNRAVARIANT metric tensor (from inverse properties)
        dgt = -gn*[Et,Ft;Ft,Gt]*gn;
        dgp = -gn*[Ep,Fp;Fp,Gp]*gn;
        
        dgt2 = -gn*[Et2,Ft2;Ft2,Gt2]*gn - 2*dgt*[Et,Ft;Ft,Gt]*gn;
        dgp2 = -gn*[Ep2,Fp2;Fp2,Gp2]*gn - 2*dgp*[Ep,Fp;Fp,Gp]*gn;
        dgtp = -gn*[Etp,Ftp;Ftp,Gtp]*gn - 2*dgt*[Ep,Fp;Fp,Gp]*gn;
            
%       Shear transverse
%       Contravariant basis vectors, raise index 
        c1 = dxt(:,i2,j)*gn(1,1) + dxp(:,i2,j)*gn(1,2);
        c2 = dxt(:,i2,j)*gn(2,1) + dxp(:,i2,j)*gn(2,2);
        
%       For derivatives of tau/q, we need dc/dtht and dc/dphi
%       Chain rule across raise index operation
        c1t = dxt2(:,i2,j)*gn(1,1) + dxtp(:,i2,j)*gn(2,1)...
            + dxt(:,i2,j)*dgt(1,1) + dxp(:,i2,j)*dgt(2,1);
        c1p = dxtp(:,i2,j)*gn(1,1) + dxp2(:,i2,j)*gn(2,1)...
            + dxt(:,i2,j)*dgp(1,1) + dxp(:,i2,j)*dgp(2,1);
        c2t = dxt2(:,i2,j)*gn(1,2) + dxtp(:,i2,j)*gn(2,2)...
            + dxt(:,i2,j)*dgt(1,2) + dxp(:,i2,j)*dgt(2,2);
        c2p = dxtp(:,i2,j)*gn(1,2) + dxp2(:,i2,j)*gn(2,2)...
            + dxt(:,i2,j)*dgp(1,2) + dxp(:,i2,j)*dgp(2,2);
        
%       Laplace-Beltrami of mean curvature. Three chain rules with 4 things to sum each
%       First: deriv of J
        LBk = 1/J(i2,j)*(Jt*gn(1,1)*kt + Jt*gn(1,2)*kp + Jp*gn(2,1)*kt + Jp*gn(2,2)*kp) ...
            + dgt(1,1)*kt + dgt(1,2)*kp + dgp(2,1)*kt + dgp(2,2)*kp ...
            + gn(1,1)*kt2 + gn(1,2)*ktp + gn(2,1)*ktp + gn(2,2)*kp2;
        
%       !!!!!!!!!! HACKY FIRST TIME STEP
        if(i == 1)
            JR(i2,j) = J(i2,j);
            c1R(:,i2,j) = c1;
            c2R(:,i2,j) = c2;
            kR(i2,j) = k;
            kdR(1,i2,j) = kt;
            kdR(2,i2,j) = kp;
            kd2R(1,i2,j)= kt2;
            kd2R(2,i2,j)= kp2;
            kd2R(3,i2,j)= ktp;
            
            c1tR(:,i2,j) = c1t;
            c1pR(:,i2,j) = c1p;
            c2tR(:,i2,j) = c2t;
            c2pR(:,i2,j) = c2p;
            
            gnR(:,:,i2,j) = gn;
            dgtR(:,:,i2,j) = dgt;
        end
        
%       Christoffels (Page 156-157 in Mollman if you want to check, but seem to be right)
        c111 = dot(dxt2(:,i2,j),c1);
        c112 = dot(dxtp(:,i2,j),c1);
        c122 = dot(dxp2(:,i2,j),c1);
        c222 = dot(dxp2(:,i2,j),c2);
        c221 = dot(dxtp(:,i2,j),c2);
        c211 = dot(dxt2(:,i2,j),c2);
        
%       Christoffel derivatives: need 5 each
        c111t = dot(dxt3(:,i2,j) ,c1) + dot(dxt2(:,i2,j),c1t);
        c112t = dot(dxt2p(:,i2,j),c1) + dot(dxtp(:,i2,j),c1t);
        c122t = dot(dxtp2(:,i2,j),c1) + dot(dxp2(:,i2,j),c1t);
        c221t = dot(dxt2p(:,i2,j),c2) + dot(dxtp(:,i2,j),c2t);
        c222t = dot(dxtp2(:,i2,j),c2) + dot(dxp2(:,i2,j),c2t);
        
        c111p = dot(dxt2p(:,i2,j),c1) + dot(dxt2(:,i2,j),c1p);
        c112p = dot(dxtp2(:,i2,j),c1) + dot(dxtp(:,i2,j),c1p);
        c211p = dot(dxt2p(:,i2,j),c2) + dot(dxt2(:,i2,j),c2p);
        c221p = dot(dxtp2(:,i2,j),c2) + dot(dxtp(:,i2,j),c2p);
        c222p = dot(dxp3(:,i2,j) ,c2) + dot(dxp2(:,i2,j),c2p);
          
%       New Helfrich bending
        fb = Eb*(2*LBk + (2*k + c0)*(2*k*k - 2*kG  - c0*k));

%%      In plane Tension
%       Deformation gradient tensor
        Fd = dxt(:,i2,j)*c1R(:,i2,j)' + dxp(:,i2,j)*c2R(:,i2,j)';
        
%       Surface projection operator
        P = eye(3) - nk(:,i2,j)*nk(:,i2,j)';
        
%       Cauchy-Green
        V2 = Fd*Fd';
        
%       Principal strains
        [ev,lams] = eig(V2);
        es = sqrt(diag(lams));
        ps1(i2,j) = es(3);
        ps2(i2,j) = es(2);
        Ar(i2,j) = es(3)*es(2);
        Sh(i2,j) = es(3)/es(2);
%       Normalized eigenvectors        
        ev1 = ev(:,1)/norm(ev(:,1));
        ev2 = ev(:,2)/norm(ev(:,2));
        
%       Strain invariants
        I1 = es(3)^2 + es(2)^2 - 2;
        I2 = es(3)^2*es(2)^2 - 1;
        
%       Covariant in plane deformation gradient tensor !!!! just a check
        Fdv = [dxt(:,i2,j)'*Fd*dxt(:,i2,j),dxt(:,i2,j)'*Fd*dxp(:,i2,j)
               dxp(:,i2,j)'*Fd*dxt(:,i2,j),dxp(:,i2,j)'*Fd*dxp(:,i2,j)];
%       Contravariant
        Fdn = [c1'*Fd*c1,c1'*Fd*c2
               c2'*Fd*c1,c2'*Fd*c2];
           
%       In plane strain tensor
        eps = (Fdn*Fdv' - eye(2))/2;
        
%       In plane tension
        tau = 0.5*(B/(es(3)*es(2))*(I1 + 1)*V2 ... % Numerical instability here: C increases I2 when it's at ~1e-16 to ~1e-14
            +  es(3)*es(2)*(C*I2-B)*P);
        
%       Principal tensions (not really needed, but could be good to check)
        taup1(i2,j) = es(3)/es(2)*(B/2*(2*I1+1)) + es(3)*es(2)*(-B/2+C/2*I2);
        taup2(i2,j) = es(2)/es(3)*(B/2*(2*I1+1)) + es(3)*es(2)*(-B/2+C/2*I2);

%       Matrix in surface coordinates (contravariant)
        tau11 = c1'*tau*c1;
        tau12 = c1'*tau*c2;
        tau21 = c2'*tau*c1;
        tau22 = c2'*tau*c2;
        
%       Now for the derivatives: first derivative of F/V via chain rule
        dFdt = dxt2(:,i2,j)*c1R(:,i2,j)'  + dxtp(:,i2,j)*c2R(:,i2,j)' ...
            + dxt(:,i2,j)*c1tR(:,i2,j)' + dxp(:,i2,j)*c2tR(:,i2,j)';
        
        dFdp = dxtp(:,i2,j)*c1R(:,i2,j)'  + dxp2(:,i2,j)*c2R(:,i2,j)' ...
            + dxt(:,i2,j)*c1pR(:,i2,j)' + dxp(:,i2,j)*c2pR(:,i2,j)';
        
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
        dI1t = 2*es(3)*des1t + 2*es(2)*des2t;
        dI1p = 2*es(3)*des1p + 2*es(2)*des2p;
        dI2t = (2*es(3)*es(2)^2)*des1t + (es(3)^2*2*es(2))*des2t;
        dI2p = (2*es(3)*es(2)^2)*des1p + (es(3)^2*2*es(2))*des2p;
        
%       Derivatives of projection tensor (double neg in 2nd part)
        dPt = -dnt*-nk(:,i2,j)' + nk(:,i2,j)*dnt';
        dPp = -dnp*-nk(:,i2,j)' + nk(:,i2,j)*dnp';
        
%       And finally the big one: derivatives of tau, separated via chain
%       rule, starting with es(3), then es(2)
        dtaut = (-B/(2*es(3)^2*es(2))*(I1+1)*V2 + es(2)/2*(C*I2-B)*P)*des1t ...
              + (-B/(2*es(3)*es(2)^2)*(I1+1)*V2 + es(3)/2*(C*I2-B)*P)*des2t ...
              + 0.5*B/(es(3)*es(2))*V2*dI1t ... % Invariants
              + 0.5*es(3)*es(2)*C*P*dI2t ...
              + 0.5*B/(es(3)*es(2))*(I1+1)*dV2t ...% Tensors
              + 0.5*es(3)*es(2)*(I2-B)*dPt;
          
        dtaup = (-B/(2*es(3)^2*es(2))*(I1+1)*V2 + es(2)/2*(C*I2-B)*P)*des1p ...
              + (-B/(2*es(3)*es(2)^2)*(I1+1)*V2 + es(3)/2*(C*I2-B)*P)*des2p ...
              + 0.5*B/(es(3)*es(2))*V2*dI1p ... % Invariants
              + 0.5*es(3)*es(2)*C*P*dI2p ...
              + 0.5*B/(es(3)*es(2))*(I1+1)*dV2p ...% Tensors
              + 0.5*es(3)*es(2)*(I2-B)*dPp;
          
%       Now put into dtau matrix
        dtauab = [c1t'*tau*c1 + c1'*dtaut*c1 + c1'*tau*c1t, c1t'*tau*c2 + c1'*dtaut*c2 + c1'*tau*c2t 
                  c2p'*tau*c1 + c2'*dtaup*c1 + c2'*tau*c1p, c2p'*tau*c2 + c2'*dtaup*c2 + c2'*tau*c2p];
        
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
%         cq = dq(1) + dq(2) + c111*q1 + c112*q2 + c221*q1 + c222*q2;
        
%       And finally the equilibrium equations to get the forces
        fab(1,i2,j) = -cvt; % bm(1,1)*q1 + bm(2,1)*q2 - cvt;
        fab(2,i2,j) = -cvp; % bm(1,2)*q1 + bm(2,2)*q2 - cvp;
        fab(3,i2,j) = -tau11*bv(1,1) - tau12*bv(1,2) - tau21*bv(2,1) - tau22*bv(2,2)...
                   + fb;% - cq;
        
%       These are in terms of (contravariant) surface vectors, put them into Cartesian
        myf(:,i2,j) = fab(1,i2,j)*dxt(:,i2,j) + fab(2,i2,j)*dxp(:,i2,j) + fab(3,i2,j)*-nk(:,i2,j);
        
    end
end

% Some imaginary components may slip in, so let's make them real
myf = real(myf);

%% Post processing stuff
%%% Note: the important things to track here are strains, which are, for
%%% area, J./JR and for shear SH
%   Interpolate these to physical positions
    x1 = real(SpHReconst(x1c,Yr));
    x2 = real(SpHReconst(x2c,Yr));
    x3 = real(SpHReconst(x3c,Yr));
    
%   Spectra
    for n=0:p
        if(n <= p)
            Ex(n+1,(i-1)/incr + 1) = norm(norm([x1c(n^2+1:(n+1)^2),x2c(n^2+1:(n+1)^2),x3c(n^2+1:(n+1)^2)]));
%           Ex(n+1,(i-1)/incr + 1) = norm(x3c(n^2+1:(n+1)^2));
        end
    end
    
    Jmax((i-1)/incr + 1) = max(max(J./JR));
    Jmin((i-1)/incr + 1) = min(min(J./JR));
    pst1((i-1)/incr + 1) = max(max(ps1));
    pst2((i-1)/incr + 1) = min(min(ps2));
%   Sh = lam1/lam2, which is a measure of shear strain according to Peng
    Shmax((i-1)/incr + 1) = max(max(Sh));
    Shmin((i-1)/incr + 1) = min(min(Sh));
    
    clf;
    sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])
%%%%%%%%% UNCOMMENT BELOW
%   Plot this timestep

    surf(squeeze(xf(1,:,:)),squeeze(xf(2,:,:)),squeeze(xf(3,:,:)), ...
        Sh,'edgecolor','none', ...
     'FaceAlpha',0.95,'FaceLighting','gouraud')
    %J./JR - 1,'edgecolor','none', ... % NOTE THAT J/JR = es(3)*es(2)!!
    %ps2,'edgecolor','none', ...
    set(gca,'nextplot','replacechildren','visible','off')
    colorbar
    shading interp;
%     caxis([1,1.6]);
    % Top down
    % Side
%     view(0,0);
%     axis([-2,2,0,2,-2,2])
%     pbaspect([1,.5,1])
    view(45,45);
    axis([-2+real(x1c(1)/(2*sqrt(pi))),2+real(x1c(1)/(2*sqrt(pi))),-2+real(x2c(1)/(2*sqrt(pi))),2+real(x2c(1)/(2*sqrt(pi))),-2+real(x3c(1)/(2*sqrt(pi))),2+real(x3c(1)/(2*sqrt(pi)))])
    pbaspect([1,1,1])

%%%%%%%%% UNCOMMENT ABOVE

% [cent,rad, angs]=ellipsoid_fit_new([reshape(x1,[10201,1]),reshape(x2,[10201,1]),reshape(x3,[101*101,1])]);
% % Dij((i-1)/incr + 1) = (rad(1)-rad(3))/(rad(1) + rad(3));
% 
% elx = vertcat(x1(:,1),flip(x1(:,51)));
% elz = vertcat(x3(:,1),flip(x3(:,51)));
% elxa((i-1)/incr + 1,:) = elx;
% elza((i-1)/incr + 1,:) = elz;
% rs = sqrt(elx.^2 + elz.^2);
% Dij((i-1)/incr + 1) = (max(rs)-min(rs))/(max(rs) + min(rs));
% incl((i-1)/incr + 1) = atan2(abs(angs(3,1)),abs(angs(1,1)))/4;
% incl(1) = 1/4;
% clf
% plot(thtf, taup1(:,1))
% plot(thtf,J(:,1)./JR(:,1)- J2(:,1)./JR(:,1))
% axis([0,pi,-5e-4,5e-4])


% clf
% %     surf(squeeze(xf(1,:,:)),squeeze(xf(2,:,:)),squeeze(xf(3,:,:)), ...
% %         Sh,'edgecolor','none', ...
% %      'FaceAlpha',1,'FaceLighting','gouraud')
%     surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
%          'FaceAlpha',0.75,'FaceLighting','gouraud')
%     lightangle(gca,150,50)
%     set(gca,'nextplot','replacechildren','visible','off')
%     caxis([.75, 1.5]);
%     view(90,45);
%     axis([-2,2,-2,2,-2,2])
%     pbaspect([1,1,1])
%     set(gcf, 'color', 'white');

    drawnow

% disp(max(max(max(abs(myf))))*2/B);
% % Capture the plot as an image 
% h = figure(1);
% frame = getframe(h); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(frame.cdata,256,'nodither');
% % Write to the GIF File 
% if i == 1
%   filename = 'gif2.gif';
%   imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
% else 
%   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
% end
end

% Reads the txt file output from the fortran code

% Read in raw data
% fID = fopen('fortran/x25W_2.txt');
% fID = fopen('fortran/dat/cm2v98eb03runs/xp16_1_005.txt');
% fID = fopen('fortran/x_0_03_997.txt');
fID = fopen('pap_dat/MeshIndv97/x14Ca1TT.txt');
a = fscanf(fID,'%f');
fclose(fID);
    
% fID = fopen('fortran/u_x25W_2.txt');
% fID = fopen('fortran/dat/cm2v98eb03runs/u_xp16_1_005.txt');
% fID = fopen('fortran/f_0_03_997.txt');
fID = fopen('pap_dat/MeshIndv97/u_x14Ca1TT.txt');
us = fscanf(fID,'%f');
fclose(fID);

% % fID = fopen('fortran/f_x25W_2.txt');
% fID = fopen('fortran/dat/f_x8Ca1TT.txt');
% % fID = fopen('pap_dat/TaylorValidation/Lam1/f_x005.txt');
% fs = fscanf(fID,'%f');
% fclose(fID);

% Get info about a  

% Total time steps, including initialization
tts = a(end) +                                 0;
incr = 50;
xtop = zeros(floor(tts/incr),1);
ytop = zeros(floor(tts/incr),1);
ztop = zeros(floor(tts/incr),1);
Dij = ytop;
incl = Dij;

% Number of values without time steps
lent = length(a) - tts;

% Number total values in a time step
lent1 = lent/(tts);
% length of 1 collection in one time step
tot = lent1/6;

% Order of spherical harmonics
p = sqrt(tot) - 1;
Ex = zeros(p+1,floor(tts/incr));
Eu = zeros(p+1,floor(tts/incr));

% Evaluation of spherical harmonics for interpolation
tmpt = linspace(0,pi,100);tmpp = linspace(0,2*pi,101);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmTNew(p,tmpth,tmpph);

tq = linspace(0,pi,15);pq = tq;
[ttq,ppq] = meshgrid(tq,pq);
Yqv = SpHarmTNew(p,ttq,ppq);
% Spherical harmonic evaluated at right hand side of sphere
% Ytrc = SpHarmTNew(p,pi/2,0);
Ytrc = SpHarmTNew(p,0,0);
% Time
ts = .005;
t = zeros(floor(tts/incr),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aa = size(fs);
% aa = aa(1)/tts;
% nttf = sqrt(aa/6);
% % nttf = sqrt(aa/2);
% nppf = 2*nttf;
% fmns = zeros(3,nttf,nppf,floor(tts/incr));
% % fmns = zeros(nttf,nppf,floor(tts/incr));
% 
% [xs,wg] = lgwt(nttf,-1,1);
% tht = acos(xs);
% dphi = 2*pi/nppf;
% phi = 0:dphi:dphi*(nppf-1)';
% 
% [ph,th] = meshgrid(phi,tht);
% Ytfs = SpHarmTNew(p,th,ph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Actual plotting
disp('Start!')
% Do timesteps
for i = 1:incr:tts
    t((i-1)/incr + 1) = i*ts - ts;
    
%   Column vector of all data in timestep
    raw = a((i-1)*lent1  + i:i*lent1 + i-1);
    
%   Individual directional coefficients, not distinct between real/imag)
    x1t = raw(1:3:end);
    x2t = raw(2:3:end);
    x3t = raw(3:3:end);
    
%   Real and imaginary separation
    x1c = x1t(1:tot) + x1t(tot+1:end)*1i;
    x2c = x2t(1:tot) + x2t(tot+1:end)*1i;
    x3c = x3t(1:tot) + x3t(tot+1:end)*1i;
    
%     x1c((p-1)*(p-1)+1:p*p) = 0;
%     x2c((p-1)*(p-1)+1:p*p) = 0;
%     x3c((p-1)*(p-1)+1:p*p) = 0;
    
%   Interpolate these to physical positions
    x1 = real(SpHReconst(x1c,Yr));
    x2 = real(SpHReconst(x2c,Yr));
    x3 = real(SpHReconst(x3c,Yr));
%   Same procedure for velocities
%   Column vector of all data in timestep
    raw = us((i-1)*lent1  + i:i*lent1 + i-1);
    
%   Individual directional coefficients, not distinct between real/imag)
    u1t = raw(1:3:end);
    u2t = raw(2:3:end);
    u3t = raw(3:3:end);
    
%   Real and imaginary separation
    u1c = u1t(1:tot) + u1t(tot+1:end)*1i;
    u2c = u2t(1:tot) + u2t(tot+1:end)*1i;
    u3c = u3t(1:tot) + u3t(tot+1:end)*1i;
    for n=0:p
        if(n <= p)
            Ex(n+1,(i-1)/incr + 1) = norm(norm([x1c(n^2+1:(n+1)^2),x2c(n^2+1:(n+1)^2),x3c(n^2+1:(n+1)^2)]));
%           Ex(n+1,(i-1)/incr + 1) = norm(x3c(n^2+1:(n+1)^2));
        end
%         Eu(n+1,(i-1)/incr + 1) = norm(norm([u1c(n^2+1:(n+1)^2),u2c(n^2+1:(n+1)^2),u3c(n^2+1:(n+1)^2)]));
        Eu(n+1,(i-1)/incr + 1) = norm(u1c(n^2+1:(n+1)^2));
    end
    
%   Interpolate these to physical positions
    u1 = real(SpHReconst(u1c,Yqv));
    u2 = real(SpHReconst(u2c,Yqv));
    u3 = real(SpHReconst(u3c,Yqv));
    xq1 = real(SpHReconst(x1c,Yqv,p));
    xq2 = real(SpHReconst(x2c,Yqv,p));
    xq3 = real(SpHReconst(x3c,Yqv,p));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    fmns(:,:,(i-1)/incr + 1) = reshape(us((i-1)*aa  + 1:i*aa), nttf, nppf);
%     fmns(:,:,:,(i-1)/incr + 1) = reshape(fs((i-1)*aa  + 1:i*aa), 3, nttf, nppf);
%     clf
%     plot(fmns(:,1,(i-1)/incr + 1))
%     plot(fmns(:,1,(i-1)/incr + 1)- fmns(:,1,1))
%     axis([0,18,0,3e-2])
%     hold on
%     plot(fmns(:,1,1))
%     xf1 = real(SpHReconst(x1c,Ytfs,p));
%     xf2 = real(SpHReconst(x2c,Ytfs,p));
%     xf3 = real(SpHReconst(x3c,Ytfs,p));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clf;
    sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])
%   Plot this timestep
    h1 = subplot(2,1,1);
    set(h1, 'Units', 'normalized');
    set(h1, 'Position', [-.1, 0.5, 1.15, .6]);

    
%     surf(xf1,xf2,xf3,squeeze(fmns(1,:,:,(i-1)/incr + 1)./fmns(1,:,:,2)),'edgecolor','none')
    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
         'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    set(gca,'nextplot','replacechildren','visible','off')
    % Top down
%     view(0,90);
    % Side
    view(0,0);
    axis([-2,2,0,2,-2,2])
    pbaspect([1,.5,1])
    hold on
    xtop((i-1)/incr + 1) = real(SpHReconst(x1c,Ytrc));
    ytop((i-1)/incr + 1) = real(SpHReconst(x2c,Ytrc)); 
    ztop((i-1)/incr + 1) = real(SpHReconst(x3c,Ytrc)); 
    scatter3(xtop((i-1)/incr + 1),ytop((i-1)/incr + 1),ztop((i-1)/incr + 1),75,'go','filled');
    quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1*2,[1,numel(xq1)]),reshape(u2*2,[1,numel(xq1)]),reshape(u3*2,[1,numel(xq1)]),'b','AutoScale','off')

    h2 = subplot(2,1,2);
    set(h2, 'Units', 'normalized');
    set(h2, 'Position', [0.05, 0, 1, .7]);

    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
          'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    
    set(gca,'nextplot','replacechildren','visible','off')
    set(gcf, 'color', 'white');
    view(45,45);
%     view(0,90);
    axis([-2,2,-2,2,-2,2])
    pbaspect([1,1,1])
    hold on
    scatter3(real(SpHReconst(x1c,Ytrc)),real(SpHReconst(x2c,Ytrc)),real(SpHReconst(x3c,Ytrc)),75,'go','filled');
    quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1,[1,numel(xq1)]),reshape(u2,[1,numel(xq1)]),reshape(u3,[1,numel(xq1)]),'b')

%     clf
% %     loglog(Ex8(:,(i-1)/incr +  1),'o')
% %     hold on
% %     loglog(Ex10(:,(i-1)/incr + 1),'x')
% %     loglog(Ex12(:,(i-1)/incr + 1),'p')
% %     loglog(Ex14(:,(i-1)/incr + 1),'s')
% %     loglog(Ex16(:,(i-1)/incr + 1),'.')
%     semilogy(Ex(:,(i-1)/incr + 1),'^')
% %     axis([2,q,1e-10,1e-1])
% %     ytop((i-1)/incr + 1) = u1(8,1);

[cent,rad, angs]=ellipsoid_fit_new([reshape(x1,[10100,1]),reshape(x2,[10100,1]),reshape(x3,[100*101,1])]);
% Dij((i-1)/incr + 1) = (rad(1)-rad(3))/(rad(1) + rad(3));

elx = vertcat(x1(:,1),x1(:,51));
elz = vertcat(x3(:,1),x3(:,51));
% elxa((i-1)/incr + 1,:) = elx;
% elza((i-1)/incr + 1,:) = elz;
rs = sqrt(elx.^2 + elz.^2);
Dij((i-1)/incr + 1) = (max(rs)-min(rs))/(max(rs) + min(rs));
incl((i-1)/incr + 1) = atan2(abs(angs(3,1)),abs(angs(1,1)))/4;
incl(1) = 1/4;

    drawnow
% Capture the plot as an image 
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

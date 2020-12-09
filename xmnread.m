% Reads the txt file output from the fortran code

% Read in raw data
fID = fopen('fortran/xmn1.txt');
a = fscanf(fID,'%f');
fclose(fID);

fID = fopen('fortran/u_xmn1.txt');
us = fscanf(fID,'%f');
fclose(fID);

% Get info about a

% Total time steps, including initialization
tts = a(end) +                                0;
incr = 20;
xtop = zeros(floor(tts/incr),1);
ytop = zeros(floor(tts/incr),1);
ztop = zeros(floor(tts/incr),1);

% Number of values without time steps
lent = length(a) - tts;

% Number total values in a time step
lent1 = lent/(tts);

% length of 1 collection in one time step
tot = lent1/6;

% Order of spherical harmonics
p = sqrt(tot) - 1;
Ex = zeros(p+1,floor(tts/incr));

% Evaluation of spherical harmonics for interpolation
tmpt = linspace(0,pi,100);tmpp = linspace(0,pi,100);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmTNew(p,tmpth,tmpph);

tq = linspace(0,pi,15);pq = tq;
[ttq,ppq] = meshgrid(tq,pq);
Yqv = SpHarmTNew(p,ttq,ppq);
% Spherical harmonic evaluated at right hand side of sphere
% Ytrc = SpHarmTNew(p,pi/2,0);
Ytrc = SpHarmTNew(p,0,0);

disp('Start!')
% Do timesteps
for i = 1:incr:tts
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
        Ex(n+1,(i-1)/incr + 1) = norm(norm([x1c(n^2+1:(n+1)^2),x2c(n^2+1:(n+1)^2),x3c(n^2+1:(n+1)^2)]));
%         Ex(n+1,(i-1)/incr + 1) = norm(x1c(n^2+1:(n+1)^2));
    end
    
%   Interpolate these to physical positions
    u1 = real(SpHReconst(u1c,Yqv));
    u2 = real(SpHReconst(u2c,Yqv));
    u3 = real(SpHReconst(u3c,Yqv));
    xq1 = real(SpHReconst(x1c,Yqv));
    xq2 = real(SpHReconst(x2c,Yqv));
    xq3 = real(SpHReconst(x3c,Yqv));
    
    clf;
%   Plot this timestep
    h1 = subplot(2,1,1);
    set(h1, 'Units', 'normalized');
    set(h1, 'Position', [-.1, 0.5, 1.15, .6]);

    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
         'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    set(gca,'nextplot','replacechildren','visible','off')
    % Top down
    % view(0,90);
    % Side
    view(0,0);
    axis([-2,2,0,2,-2,2])
    pbaspect([1,.5,1])
    hold on
    xtop((i-1)/incr + 1) = real(SpHReconst(x1c,Ytrc));
    ytop((i-1)/incr + 1) = real(SpHReconst(x2c,Ytrc)); 
    ztop((i-1)/incr + 1) = real(SpHReconst(x3c,Ytrc)); 
    scatter3(xtop((i-1)/incr + 1),ytop((i-1)/incr + 1),ztop((i-1)/incr + 1),75,'go','filled');
    quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1,[1,numel(xq1)]),reshape(u2,[1,numel(xq1)]),reshape(u3,[1,numel(xq1)]),'b')

    h2 = subplot(2,1,2);
    set(h2, 'Units', 'normalized');
    set(h2, 'Position', [0.05, 0, 1, .7]);

    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ... %!!!!!!!!!!!!! make full cell?
          'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    set(gca,'nextplot','replacechildren','visible','off')
    set(gcf, 'color', 'white');
    view(45,45);
    axis([-2,2,0,2,-2,2])
    pbaspect([1,.5,1])
    hold on
    scatter3(real(SpHReconst(x1c,Ytrc)),real(SpHReconst(x2c,Ytrc)),real(SpHReconst(x3c,Ytrc)),75,'go','filled');

    
%     clf
% %     loglog(Ex8(:,(i-1)/incr + 1),'o')
% %     hold on
% %     loglog(Ex10(:,(i-1)/incr + 1),'x')
% %     loglog(Ex12(:,(i-1)/incr + 1),'p')
% %     loglog(Ex14(:,(i-1)/incr + 1),'s')
% %     loglog(Ex16(:,(i-1)/incr + 1),'.')
%     loglog(Ex(:,(i-1)/incr + 1),'^')
% %     axis([2,14,1e-8,1e1])
%     xtop((i-1)/incr + 1) = u1(1,1);
    drawnow
end
% Reads the txt file output from the fortran code

% Read in raw data
fID = fopen('fortran/xmn.txt');
a = fscanf(fID,'%f');
fclose(fID);

% Get info about a

% Total time steps, including initialization
tts = a(end) + 1;

% Number of values without time steps
lent = length(a) - tts;

% Number total values in a time step
lent1 = lent/tts;

% length of 1 collection in one time step
tot = lent1/6;

% Order of spherical harmonics
p = sqrt(tot) - 1;

% Evaluation of spherical harmonics for interpolation
tmpt = linspace(0,pi,100);tmpp = linspace(0,pi,100);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmTNew(p,tmpth,tmpph);
% Spherical harmonic evaluated at right hand side of sphere
% Ytrc = SpHarmTNew(p,pi/2,0);
Ytrc = SpHarmTNew(p,0,0);

disp('Start!')
% Do timesteps
for i = 1:15:tts
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
    scatter3(real(SpHReconst(x1c,Ytrc)),real(SpHReconst(x2c,Ytrc)),real(SpHReconst(x3c,Ytrc)),75,'go','filled');

    % Just from stress - magenta
    % quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua(1,:,:),[1,numel(x(1,:,:))]),reshape(ua(2,:,:),[1,numel(x(1,:,:))]),reshape(ua(3,:,:),[1,numel(x(1,:,:))]),'m')
    % Force - green
    % quiver3(reshape(xf(1,:,:),[1,numel(xf(1,:,:))]),reshape(xf(2,:,:),[1,numel(xf(1,:,:))]),reshape(xf(3,:,:),[1,numel(xf(1,:,:))]),reshape(real(myf(1,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(2,:,:)),[1,numel(xf(1,:,:))]),reshape(real(myf(3,:,:)),[1,numel(xf(1,:,:))]),'g')
    % Stress and fluid - blue
    % quiver3(reshape(x(1,:,:),[1,numel(x(1,:,:))]),reshape(x(2,:,:),[1,numel(x(1,:,:))]),reshape(x(3,:,:),[1,numel(x(1,:,:))]),reshape(ua11(1,:,:),[1,numel(x(1,:,:))]),reshape(ua11(2,:,:),[1,numel(x(1,:,:))]),reshape(ua11(3,:,:),[1,numel(x(1,:,:))]),'b')

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
    drawnow
end
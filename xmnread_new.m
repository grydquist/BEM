% Reads the txt file output from the fortran code
fclose all;
% Get total timesteps outputted

dir = 'fortran/dat/TWZp16F340pN/';
fid = fopen(strcat(dir,'maxdt'));
tts = str2double(fgetl(fid));
fclose(fid);
tts = floor(tts);

% Read Param file to get simulation info
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
    end
    tline = fgetl(fid);
end
fclose(fid);

% How many timesteps to skip
incr = 25;
% Round down to fit w/ dt_inc
incr = incr - mod(incr,dt_inc);
if(incr == 0); incr = dt_inc; end
% incr = dt_inc*10;




tsteps = floor(tts/incr) + 1;
xtop = zeros(tsteps,1);
ytop = zeros(tsteps,1);
ztop = zeros(tsteps,1);
Dij = ytop;
incl = Dij;

% Number total values in a time step  
lent1 = 6*(p+1)^2;

% Order of spherical harmonics
tot = (p+1)^2;
Ex = zeros(p+1,tsteps);
Eu = zeros(p+1,tsteps);

% Evaluation of spherical harmonics for interpolation
tmpt = linspace(0,pi,101);tmpp = linspace(0,2*pi,101);
[tmpph,tmpth] = meshgrid(tmpp,tmpt);
Yr = SpHarmTNew(p,tmpth,tmpph);

tq = linspace(0,pi,15);pq = tq;
[ttq,ppq] = meshgrid(tq,pq);
Yqv = SpHarmTNew(p,ttq,ppq);
% Spherical harmonic evaluated at right hand side of sphere
% Ytrc = SpHarmTNew(p,pi/2,0);
Ytrc = SpHarmTNew(p,pi/2,0);
% Time
t = zeros(tsteps,1);

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
    
%     x1c((p-1)*(p-1)+1:p*p) = 0;
%     x2c((p-1)*(p-1)+1:p*p) = 0;
%     x3c((p-1)*(p-1)+1:p*p) = 0;
    
%   Interpolate these to physical positions
    x1 = real(SpHReconst(x1c,Yr));
    x2 = real(SpHReconst(x2c,Yr));
    x3 = real(SpHReconst(x3c,Yr));
%   Same procedure for velocities
%   Read in data file of current timestep, making exception for the first
    if(i == 1)
        file = 'u_00001';
    else
        file = 'u_';
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

    clf;
%   Plot this timestep
    h1 = subplot(2,1,1);
    sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])
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
%     semilogy(Eu(:,(i-1)/incr + 1),'^')
% %     axis([2,q,1e-10,1e-1])
% %     ytop((i-1)/incr + 1) = u1(8,1);

[cent,rad, angs]=ellipsoid_fit_new([reshape(x1,[10201,1]),reshape(x2,[10201,1]),reshape(x3,[101*101,1])]);
% Dij((i-1)/incr + 1) = (rad(1)-rad(3))/(rad(1) + rad(3));

elx = vertcat(x1(:,1),x1(:,51));
elz = vertcat(x3(:,1),x3(:,51));
elxa((i-1)/incr + 1,:) = elx;
elza((i-1)/incr + 1,:) = elz;
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

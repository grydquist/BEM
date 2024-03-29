% Reads the txt file output from the fortran code
fclose all;

% Path to data directory
% dir = 'fortran/dat/TEST1/';
% dir = 'fortran/dat/new8bnd25/new8bnd252/';
% dir = 'fortran/dat/new16bnd25/new16bnd252/';
dir = 'fortran/dat/TEST/cell_1/';
dir = 'fortran/dat/P13S10H45/cell_8/';
% dir = 'fortran/dat/reduces/Reduce18/';
% dir = 'fortran/dat/repa12_8_25/cell_2/';
% dir = 'fortran/dat/Param12_test/cell_2/';
% dir = 'fortran/dat/Hct20_12_8/cell_6/';
% 
% dir = 'fortran/dat/SR10Hct25/cell_3/';
% dir = 'fortran/dat/p16Ca02/cell_1/';
% dir = 'fortran/dat/NA500H45/cell_1/';
% dir = 'fortran/dat/IndNewSR25/cell_1/';
% dir = 'fortran/dat/p12S100H45/cell_2/';

% dir = 'fortran/dat/ERAZURE/cell_1/';
% % dir = 'fortran/dat/bend_only/bend_only2/';
% dir = 'fortran/dat/rbcsxtn/rbcsxtn4/';

% Max timestep
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

% Time increment to plot results
tincr = .1;

incr = floor(tincr/ts);%500;
% Round down to fit w/ dt_inc
incr = incr - mod(incr,dt_inc);
if(incr == 0); incr = dt_inc; end
 
% Tracking point along top
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


gf=1;
ntf = gf*(p+1);
npf = gf*2*(p+1);
dphif = 2*pi/npf;
pq = 0:dphif:dphif*(npf-1)';
[xsf,wgf] = lgwt(ntf,-1,1);
tq = acos(xsf);
[ttq,ppq] = meshgrid(tq,pq);

Yqv = SpHarmTNew(p,ttq,ppq);
% Spherical harmonic evaluated at top of sphere
Ytrc = SpHarmTNew(p,0,0);

% Time
t = zeros(tsteps,1);

%% Actual plotting
disp('Start!')
% Do timesteps
for i = 1:incr:tts + 1
% Current time
    t((i-1)/incr + 1) = i*ts - ts;
    
%%  Reading in the files and separating into useful quantities
%   Read in data file of current timestep, making exception for the first
    if(i == 1)
        file = 'x_00000'; %%!! inconsistent with some older files (change to 1)
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
    
%   Interpolate these to physical positions
    x1 = real(SpHReconst(x1c,Yr));
    x2 = real(SpHReconst(x2c,Yr));
    x3 = real(SpHReconst(x3c,Yr));
    
    
%   Exact same procedure for velocity
%   Read in data file of current timestep, making exception for the first
    if(i == 1)
        file = 'u_00000'; %%!! inconsistent with some older files (change to 1)
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
    
%   Interpolate these to physical positions, as well as
    u1 = real(SpHReconst(u1c,Yqv));
    u2 = real(SpHReconst(u2c,Yqv));
    u3 = real(SpHReconst(u3c,Yqv));
    xq1 = real(SpHReconst(x1c,Yqv,p));
    xq2 = real(SpHReconst(x2c,Yqv,p));
    xq3 = real(SpHReconst(x3c,Yqv,p));
    
% %   Uncomment for forces
% %   Read in data file of current timestep, making exception for the first
%     if(i == 1)
%         file = 'f_00000'; %%!! inconsistent with some older files (change to 1)
%     else
%         file = 'f_';
%         for fl = 1:(4-floor(log10(i-1)))
%             file = strcat(file,'0');
%         end
%         file = strcat(file,num2str(i-1));
%     end
%     
% %   All the data in the ts
%     fID = fopen(strcat(dir,file));
%     raw = fscanf(fID,'%f');
%     fclose(fID);
%     
% %   Individual directional coefficients, not distinct between real/imag)
%     f1t = raw(1:3:end);
%     f2t = raw(2:3:end);
%     f3t = raw(3:3:end);
%     
% %   Real and imaginary separation
%     f1c = f1t(1:tot) + f1t(tot+1:end)*1i;
%     f2c = f2t(1:tot) + f2t(tot+1:end)*1i;
%     f3c = f3t(1:tot) + f3t(tot+1:end)*1i;
% %   Interpolate these to physical positions
%     f1 = real(SpHReconst(f1c,Yqv,p));
%     f2 = real(SpHReconst(f2c,Yqv,p));
%     f3 = real(SpHReconst(f3c,Yqv,p));
    
    
%   Uncomment for spectra
%     for n=0:p
%         if(n <= p)
%             Ex(n+1,(i-1)/incr + 1) = norm(norm([x1c(n^2+1:(n+1)^2),x2c(n^2+1:(n+1)^2),x3c(n^2+1:(n+1)^2)]));
% %           Ex(n+1,(i-1)/incr + 1) = norm(x3c(n^2+1:(n+1)^2));
%         end
% %         Eu(n+1,(i-1)/incr + 1) = norm(norm([u1c(n^2+1:(n+1)^2),u2c(n^2+1:(n+1)^2),u3c(n^2+1:(n+1)^2)]));
%         Eu(n+1,(i-1)/incr + 1) = norm(u1c(n^2+1:(n+1)^2));
%     end
    

%%  Actual plotting of above quantities
    clf;
%   Plot this timestep

%   First subplot
    h1 = subplot(2,1,1);
    sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])
    set(h1, 'Units', 'normalized');
    set(h1, 'Position', [-.1, 0.5, 1.15, .6]);

%   Surface    
    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
         'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,25)
    set(gca,'nextplot','replacechildren','visible','off')
    % Top down
%     view(0,90);
    % Side
    view(0,0);
    axis([0,10,0,5,2,5])
    pbaspect([2,1,3/5])
    hold on
    xtop((i-1)/incr + 1) = real(SpHReconst(x1c,Ytrc));
    ytop((i-1)/incr + 1) = real(SpHReconst(x2c,Ytrc)); 
    ztop((i-1)/incr + 1) = real(SpHReconst(x3c,Ytrc)); 
    
%   Plot material point
    scatter3(xtop((i-1)/incr + 1),ytop((i-1)/incr + 1),ztop((i-1)/incr + 1),75,'go','filled');
%   Plot velocity vectors
%     quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1*2,[1,numel(xq1)]),reshape(u2*2,[1,numel(xq1)]),reshape(u3*2,[1,numel(xq1)]),'b')%,'AutoScale','off')
    scatter3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]))

%   Second subplot, just a different view.
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
    axis([0,5,0,5,2,5])
    pbaspect([1,1,3/5])
    hold on
    scatter3(real(SpHReconst(x1c,Ytrc)),real(SpHReconst(x2c,Ytrc)),real(SpHReconst(x3c,Ytrc)),75,'go','filled');
    quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1,[1,numel(xq1)]),reshape(u2,[1,numel(xq1)]),reshape(u3,[1,numel(xq1)]),'b')

% Extra stuff

% [cent,rad, angs]=ellipsoid_fit_new([reshape(x1,[10201,1]),reshape(x2,[10201,1]),reshape(x3,[101*101,1])]);
% % Dij((i-1)/incr + 1) = (rad(1)-rad(3))/(rad(1) + rad(3));
% 
% elx = vertcat(x1(:,1),flip(x1(:,51)));
% elz = vertcat(x3(:,1),flip(x3(:,51)));
% % elxa((i-1)/incr + 1,:) = elx;
% % elza((i-1)/incr + 1,:) = elz;
% rs = sqrt(elx.^2 + elz.^2);
% Dij((i-1)/incr + 1) = (max(rs)-min(rs))/(max(rs) + min(rs));
% incl((i-1)/incr + 1) = atan2(abs(angs(3,1)),abs(angs(1,1)))/4;
% incl(1) = 1/4; 


clf
    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
         'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    set(gca,'nextplot','replacechildren','visible','off')
    % Top down
%     view(0,90);
    % Side
    view(0,0);
    axis([-2,2,-2,2,-2,2])
    pbaspect([1,1,1])
    
    set(gca,'nextplot','replacechildren','visible','off')
    set(gcf, 'color', 'white');
%     view(45,45);
% clf
% surf(xq1,xq2,xq3);%, f1)
% xm1=mean(mean(xq1));
% xm2=mean(mean(xq2));
% xm3=mean(mean(xq3));
% axis([xm1-1.5,xm1+1.5, xm2-1.5,xm2+1.5,xm3-1.5,xm3+1.5])
% pbaspect([1,1,1])
% view(45,45)
% title(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])

    drawnow
    if(t((i-1)/incr + 1)>.24)
    lll=1;
    end

% % Capture the plot as an image 
% h = figure(1);
% frame = getframe(h); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(frame.cdata,256,'nodither');
% % Write to the GIF File 
% if i == 1
%   filename = 'gif3.gif';
%   imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
% else 
%   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
% end
end

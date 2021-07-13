% Reads the txt file output from the fortran code
fclose all;
% Get total timesteps outputted

% dir = 'pap_dat/MeshIndNew/TT18/';
% dir = 'pap_dat/TurbRes/Ca3/HITCa3_11/';
% dir = 'fortran/dat/tmptmp';
% dir = 'fortran/dat/tmptmptmp';
dir = 'fortran/dat/tmp';
% dir = 'fortran/dat/ERAZURE7/';
% dir = 'pap_dat/TurbRes/Ca1/HITCa1_2/';
% dir = 'pap_dat/TWZ/cmplTWZ/TWZp16F90pN/';
fid = fopen(strcat(dir,'1/','maxdt'));
tts = str2double(fgetl(fid));
fclose(fid);
tts = floor(tts);

% Read Param file to get simulation info
fid = fopen(strcat(dir,'1/','Params'));
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
incr = floor(.125/ts);%500;
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
Ytrc = SpHarmTNew(p,0,0);
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
        file = 'x_00000';
    else
        file = 'x_';
        for fl = 1:(4-floor(log10(i-1)))
            file = strcat(file,'0');
        end
        file = strcat(file,num2str(i-1));
    end
    
    %% Do this for two cells!
    
%   All the data in the ts
    fID = fopen(strcat(dir,'1/',file));
    raw = fscanf(fID,'%f');
    fclose(fID);
    fID = fopen(strcat(dir,'2/',file));
    raw2 = fscanf(fID,'%f');
    fclose(fID);
    
%   Individual directional coefficients, not distinct between real/imag)
    x1t = raw(1:3:end);
    x2t = raw(2:3:end);
    x3t = raw(3:3:end);
    x1t2 = raw2(1:3:end);
    x2t2 = raw2(2:3:end);
    x3t2 = raw2(3:3:end);
    
%   Real and imaginary separation
    x1c = x1t(1:tot) + x1t(tot+1:end)*1i;
    x2c = x2t(1:tot) + x2t(tot+1:end)*1i;
    x3c = x3t(1:tot) + x3t(tot+1:end)*1i;
    x1c2 = x1t2(1:tot) + x1t2(tot+1:end)*1i;
    x2c2 = x2t2(1:tot) + x2t2(tot+1:end)*1i;
    x3c2 = x3t2(1:tot) + x3t2(tot+1:end)*1i;
    
%   Interpolate these to physical positions
    x1 = real(SpHReconst(x1c,Yr));
    x2 = real(SpHReconst(x2c,Yr));
    x3 = real(SpHReconst(x3c,Yr));
    x12 = real(SpHReconst(x1c2,Yr));
    x22 = real(SpHReconst(x2c2,Yr));
    x32 = real(SpHReconst(x3c2,Yr));
%   Same procedure for velocities
%   Read in data file of current timestep, making exception for the first
    if(i == 1)
        file = 'u_00000';
    else
        file = 'u_';
        for fl = 1:(4-floor(log10(i-1)))
            file = strcat(file,'0');
        end
        file = strcat(file,num2str(i-1));
    end
    
%   All the data in the ts
    fID = fopen(strcat(dir,'1/',file));
    raw = fscanf(fID,'%f');
    fclose(fID);
    fID = fopen(strcat(dir,'2/',file));
    raw2 = fscanf(fID,'%f');
    fclose(fID);
    
%   Individual directional coefficients, not distinct between real/imag)
    u1t = raw(1:3:end);
    u2t = raw(2:3:end);
    u3t = raw(3:3:end);
    u1t2 = raw2(1:3:end);
    u2t2 = raw2(2:3:end);
    u3t2 = raw2(3:3:end);
    
%   Real and imaginary separation
    u1c = u1t(1:tot) + u1t(tot+1:end)*1i;
    u2c = u2t(1:tot) + u2t(tot+1:end)*1i;
    u3c = u3t(1:tot) + u3t(tot+1:end)*1i;
    u1c2 = u1t2(1:tot) + u1t2(tot+1:end)*1i;
    u2c2 = u2t2(1:tot) + u2t2(tot+1:end)*1i;
    u3c2 = u3t2(1:tot) + u3t2(tot+1:end)*1i;
    
%   Interpolate these to physical positions
    u1 = real(SpHReconst(u1c,Yqv));
    u2 = real(SpHReconst(u2c,Yqv));
    u3 = real(SpHReconst(u3c,Yqv));
    xq1 = real(SpHReconst(x1c,Yqv,p));
    xq2 = real(SpHReconst(x2c,Yqv,p));
    xq3 = real(SpHReconst(x3c,Yqv,p));
    u12 = real(SpHReconst(u1c2,Yqv));
    u22 = real(SpHReconst(u2c2,Yqv));
    u32 = real(SpHReconst(u3c2,Yqv));
    xq12 = real(SpHReconst(x1c2,Yqv,p));
    xq22 = real(SpHReconst(x2c2,Yqv,p));
    xq32 = real(SpHReconst(x3c2,Yqv,p));

    clf;
%   Plot this timestep
    sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])

    
%     surf(xf1,xf2,xf3,squeeze(fmns(1,:,:,(i-1)/incr + 1)./fmns(1,:,:,2)),'edgecolor','none')
    surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
         'FaceAlpha',0.75,'FaceLighting','gouraud')
    lightangle(gca,150,50)
    hold on
    surf(x12,x22,x32,'edgecolor','none','FaceColor',[1 0 0], ...
     'FaceAlpha',0.75,'FaceLighting','gouraud')
    
    set(gca,'nextplot','replacechildren','visible','off')
    % Top down
%     view(0,90);
    % Side
    view(0,0);
    axis([-3,3,0,2,-2,2])
    pbaspect([6,2,4])
    hold on
    xtop((i-1)/incr + 1) = real(SpHReconst(x1c,Ytrc));
    ytop((i-1)/incr + 1) = real(SpHReconst(x2c,Ytrc)); 
    ztop((i-1)/incr + 1) = real(SpHReconst(x3c,Ytrc)); 
%     scatter3(xtop((i-1)/incr + 1),ytop((i-1)/incr + 1),ztop((i-1)/incr + 1),75,'go','filled');
    quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1*2,[1,numel(xq1)]),reshape(u2*2,[1,numel(xq1)]),reshape(u3*2,[1,numel(xq1)]),'b')%,'AutoScale','off')
    quiver3(reshape(xq12,[1,numel(xq12)]),reshape(xq22,[1,numel(xq22)]),reshape(xq32,[1,numel(xq12)]),reshape(u12*2,[1,numel(xq12)]),reshape(u22*2,[1,numel(xq12)]),reshape(u32*2,[1,numel(xq12)]),'b')%,'AutoScale','off')

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

elx = vertcat(x1(:,1),flip(x1(:,51)));
elz = vertcat(x3(:,1),flip(x3(:,51)));
% elxa((i-1)/incr + 1,:) = elx;
% elza((i-1)/incr + 1,:) = elz;
rs = sqrt(elx.^2 + elz.^2);
Dij((i-1)/incr + 1) = (max(rs)-min(rs))/(max(rs) + min(rs));
incl((i-1)/incr + 1) = atan2(abs(angs(3,1)),abs(angs(1,1)))/4;
incl(1) = 1/4;

% 
% clf
%     surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
%          'FaceAlpha',0.75,'FaceLighting','gouraud')
%     lightangle(gca,150,50)
%     set(gca,'nextplot','replacechildren','visible','off')
%     % Top down
% %     view(0,90);
%     % Side
%     view(0,0);
%     axis([-2,2,-2,2,-2,2])
%     pbaspect([1,1,1])
%     
%     set(gca,'nextplot','replacechildren','visible','off')
%     set(gcf, 'color', 'white');
%     view(45,45);
    
    drawnow
% % Capture the plot as an image 
% h = figure(1);
% frame = getframe(h); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(frame.cdata,256,'nodither');
% % Write to the GIF File 
% if i == 1
%   filename = 'pull3.gif';
%   imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
% else 
%   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
% end
end

% Reads the txt file output from the fortran code
fclose all;
% Get total timesteps outputted

% dir = 'pap_dat/MeshIndNew/TT18/';
% dir = 'pap_dat/TurbRes/Ca3/HITCa3_11/';
% dir = 'fortran/dat/tmptmp';
% dir = 'fortran/dat/tmptmptmp';
% dir = 'fortran/dat/tmp';
dir = 'fortran/dat/mult_cells';
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
    
    if(i == 1)
        file = 'x_00000';
    else
        file = 'x_';
        for fl = 1:(4-floor(log10(i-1)))
            file = strcat(file,'0');
        end
        file = strcat(file,num2str(i-1));
    end
%%  Loop through all files
    celln = 1;
    curfile = strcat(dir,num2str(celln),'/',file);
    clf
    
    while(isfile(curfile))
        
%       All the data in the ts
        fID = fopen(curfile);
        raw = fscanf(fID,'%f');
        fclose(fID);

%       Individual directional coefficients, not distinct between real/imag)
        x1t = raw(1:3:end);
        x2t = raw(2:3:end);
        x3t = raw(3:3:end);

%       Real and imaginary separation
        x1c = x1t(1:tot) + x1t(tot+1:end)*1i;
        x2c = x2t(1:tot) + x2t(tot+1:end)*1i;
        x3c = x3t(1:tot) + x3t(tot+1:end)*1i;

%       Interpolate these to physical positions
        x1 = real(SpHReconst(x1c,Yr));
        x2 = real(SpHReconst(x2c,Yr));
        x3 = real(SpHReconst(x3c,Yr));
        
%       Same procedure for velocities
%       Read in data file of current timestep, making exception for the first
        if(i == 1)
            file = 'u_00000';
        else
            file = 'u_';
            for fl = 1:(4-floor(log10(i-1)))
                file = strcat(file,'0');
            end
            file = strcat(file,num2str(i-1));
        end

%       All the data in the ts
        fID = fopen(strcat(dir,num2str(celln),'/',file));
        raw = fscanf(fID,'%f');
        fclose(fID);

%       Individual directional coefficients, not distinct between real/imag)
        u1t = raw(1:3:end);
        u2t = raw(2:3:end);
        u3t = raw(3:3:end);

%       Real and imaginary separation
        u1c = u1t(1:tot) + u1t(tot+1:end)*1i;
        u2c = u2t(1:tot) + u2t(tot+1:end)*1i;
        u3c = u3t(1:tot) + u3t(tot+1:end)*1i;

%       Interpolate these to physical positions
        u1 = real(SpHReconst(u1c,Yqv));
        u2 = real(SpHReconst(u2c,Yqv));
        u3 = real(SpHReconst(u3c,Yqv));
        xq1 = real(SpHReconst(x1c,Yqv,p));
        xq2 = real(SpHReconst(x2c,Yqv,p));
        xq3 = real(SpHReconst(x3c,Yqv,p));

%       Plot this timestep
        sgtitle(['time = ',num2str(t((i-1)/incr + 1)),',  iter = ',num2str(i)])

%       surf(xf1,xf2,xf3,squeeze(fmns(1,:,:,(i-1)/incr + 1)./fmns(1,:,:,2)),'edgecolor','none')
        surf(x1,x2,x3,'edgecolor','none','FaceColor',[1 0 0], ...
             'FaceAlpha',0.75,'FaceLighting','gouraud')
        lightangle(gca,150,50)
        hold on

        set(gca,'nextplot','replacechildren','visible','off')
        % Top down
    %     view(0,90);
        % Side
        view(0,0);
        axis([-6,6,0,2,-2,2])
        pbaspect([12,2,4])
        hold on
        quiver3(reshape(xq1,[1,numel(xq1)]),reshape(xq2,[1,numel(xq2)]),reshape(xq3,[1,numel(xq1)]),reshape(u1*2,[1,numel(xq1)]),reshape(u2*2,[1,numel(xq1)]),reshape(u3*2,[1,numel(xq1)]),'b')%,'AutoScale','off')

%       Read in data file of current timestep, making exception for the first
        celln = celln + 1;
        if(i == 1)
            file = 'x_00000';
        else
            file = 'x_';
            for fl = 1:(4-floor(log10(i-1)))
                file = strcat(file,'0');
            end
            file = strcat(file,num2str(i-1));
        end
        curfile = strcat(dir,num2str(celln),'/',file);
    end
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

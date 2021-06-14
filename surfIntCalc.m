%% Loading in
p = 12;
fali = 2;
q = p*fali;
np = 2*(q+1);
nt = q+1;

dphi = 2*pi/np;
phif = 0:dphi:dphi*(np-1)';
[trash,wg] = lgwt(nt,-1,1);
thtf = acos(trash);
[phf,thf] = meshgrid(phif,thtf);


xsx12 = dlmread('./fortran/tmperx');
xsx12 = reshape(xsx12,[np,nt])';
xsy12 = dlmread('./fortran/tmpery');
xsy12 = reshape(xsy12,[np,nt])';
xsz12 = dlmread('./fortran/tmperz');
xsz12 = reshape(xsz12,[np,nt])';
fsx12 = dlmread('./fortran/tmper1');
fsx12 = reshape(fsx12,[np,nt])';
fsy12 = dlmread('./fortran/tmper2');
fsy12 = reshape(fsy12,[np,nt])';
fsz12 = dlmread('./fortran/tmper3');
fsz12 = reshape(fsz12,[np,nt])';
J = dlmread('./fortran/tmperJ');
J = reshape(J,[np,nt])';
surf(xsx12,xsy12,xsz12,J./sin(thf))%sqrt(fsx12.^2+fsy12.^2+fsz12.^2))
axis([-1.5,1.5,-1,1,-1,1]);pbaspect([3,2,2])

nx = 20;
xq = linspace(-3,3,nx);
nz = 40;
zq = linspace(-3,3,nz);

[XQ,ZQ] = meshgrid(xq,zq);
bx = zeros(nx,nz);
by = bx;
bz = bx;

%% Calculate the integral (which corresponds to velocity) at each point

% Evaluation point
for i = 1:nx
    for j = 1:nz
        b = [0,0,0]';
%       Integral eval
        for i2 = 1:nt
            for j2 = 1:np
                r = [xsx12(i2,j2),xsy12(i2,j2),xsz12(i2,j2)] - [xq(i),0,zq(j)];
                f = [fsx12(i2,j2),fsy12(i2,j2),fsz12(i2,j2)];
                Gs = Gij(r);
                b = b + wg(i2)*dphi*Gs*f'*J(i2,j2)/sin(thtf(i2));
            end
        end
        
        bx(i,j) = b(1);
        by(i,j) = b(2);
        bz(i,j) = b(3);
    end
end

quiver(xq,zq,-bx',-bz')


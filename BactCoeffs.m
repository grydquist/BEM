function [x1mn,x2mn,x3mn, x1] = BactCoeffs(Y, tht, phi, ord) 
% Given the order of spherical harmonics and points to evaluate, returns x,
%y, and z coeffecients for RBC shape

% Cutoff angle for cap
h = 1;
thtc = atan(h);


[t,p] = size(tht);
thts = tht(:,1);

x1 = zeros(t,p);
x2 = x1;
x3 = x1;

% Rotate arond the phis.
for i = 1:t
    for j = 1:p
        ct = pi/2-thts(i);
        if ct>thtc
            th3 = pi-ct;
            zz = h/tan(ct);
            th4 = asin(sin(th3)*zz);
            th2 = pi - th3 - th4;
            rr = sqrt(zz^2 + 1 - 2*zz*cos(th2));
            xx = rr*cos(ct);
            xy = zz + xx;
            x3(i,j) = h + xx*tan(ct);
            
        elseif abs(ct)>thtc
            ct = abs(ct);
            th3 = pi-ct;
            zz = h/tan(ct);
            th4 = asin(sin(th3)*zz);
            th2 = pi - th3 - th4;
            rr = sqrt(zz^2 + 1 - 2*zz*cos(th2));
            xx = rr*cos(ct);
            xy = zz + xx;
            x3(i,j) = -(h + xx*tan(ct));
        else
            xy = 1;
            x3(i,j) = tan(ct);
        end
        x1(i,j) = xy*cos(phi(i,j));
        x2(i,j) = xy*sin(phi(i,j));
        
    end
end

if(nargin ==4)
    n = ord;
else
    n = length(Y)-1;
end
% Do the transforms individually
x1mn = SpT(Y,x1,tht,phi,n);
x2mn = SpT(Y,x2,tht,phi,n);
x3mn = SpT(Y,x3,tht,phi,n);

% xt1 = real(SpHReconst(x1mn,Y));
% xt2 = real(SpHReconst(x2mn,Y));
% xt3 = real(SpHReconst(x3mn,Y));
% 
% figure
% surf(xt1,xt2,xt3,'edgecolor','none');
% hold on
% surf(x1,x2,x3,'edgecolor','none');
% axis([-1 1 -1 1 -1 1]);
% pbaspect([1,1,1])



end
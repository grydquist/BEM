function [x1mn,x2mn,x3mn, x1] = RBCcoeffs(Y, tht, phi, ord) 
% Given the order of spherical harmonics and points to evaluate, returns x,
%y, and z coeffecients for RBC shape

% First define the reference shape. The x-points are defined at the
% following points
alph = 1.38581894;
a = 1;

[t,p] = size(tht);
thts = tht(:,1);

% x and y in a const phi plane. To be rotated around phi. Phi is an actual
% angle here, tht is not.

xy = a*alph*sin(thts);
z = a*alph/2*(.207+2.003*sin(thts).^2-1.123*sin(thts).^4).*cos(thts);

x1 = zeros(t,p);
x2 = x1;
x3 = x1;

% Rotate arond the phis.
for i = 1:p
    x1(:,i) = xy*cos(phi(1,i));
    x2(:,i) = xy*sin(phi(1,i));
    x3(:,i) = z;    
    
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
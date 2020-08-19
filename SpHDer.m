function [x,n,J,dxtt,dxpt] = SpHDer(fmn,Y,tht,phi)
% Calculate the derivatives of the spherical harmonics at the supplied
% points, tht and phi, as well as normal vector and Jacobian

% Order of geometric parameterization
pmn = sqrt(length(fmn)) - 1;
Ypcur = Y{1};
[nt, np] = size(squeeze(Ypcur));

% Actual locations
x = zeros(3,nt,np);

% Absolute distance from zero of a given point
r = zeros(nt,np);

% Derivatives of just harmonics at phi/tht
dxt = zeros(nt,np);
dxp = dxt;

% Total derivatives, in Cartesian
dxpt = zeros(3,nt,np);
dxtt = dxpt;

ih = 0;
it = 0;
% Loop over harmonics
for n = 0:pmn
    ih = ih+1;
    Ypcur = Y{ih};
    im = 0;
    for m = -n:n
        im = im + 1;
        it = it + 1;
%       Current value of harmonic coefficient        
        f = fmn(it);
%       No need to do the calculation if the current f is zero
        if (f == 0); continue; end
%       Current values of harmonic at degree and order
        Ym = squeeze(Ypcur(im,:,:));
        
%       Distance from zero
        r = r + f*Ym;
        
%       Get the derivative in theta and add in (two components)
        if(m > -n)
            dxt = dxt + ...
              f*(-0.5*sqrt((n+m)*(n-m+1)).*exp( 1i*phi).*squeeze(Ypcur(im-1,:,:)));
        end
        if (m < n)
            dxt = dxt + ...
              f*( 0.5*sqrt((n-m)*(n+m+1)).*exp(-1i*phi).*squeeze(Ypcur(im+1,:,:)));
        end
%       Get the derivative in phi and add in
        dxp = dxp + f*1i*m*Ym;
    end
end

% Cartesian coordinates
x(1,:,:) = r.*sin(tht).*cos(phi);
x(2,:,:) = r.*sin(tht).*sin(phi);
x(3,:,:) = r.*cos(tht);

% Derivatives in Cartesian
dxtt(1,:,:) = dxt.*sin(tht).*cos(phi) + cos(tht).*r.*cos(phi);
dxtt(2,:,:) = dxt.*sin(tht).*sin(phi) + cos(tht).*r.*sin(phi);
dxtt(3,:,:) = dxt.*cos(tht) - sin(tht).*r;

if (max(max(max(abs(imag(dxt)))))) > 1e-14
    error('Derivatives in theta are imaginary!')
end
dxtt = real(dxtt);

dxpt(1,:,:) = dxp.*sin(tht).*cos(phi) - sin(tht).*r.*sin(phi);
dxpt(2,:,:) = dxp.*sin(tht).*sin(phi) + sin(tht).*r.*cos(phi);
dxpt(3,:,:) = dxp.*cos(tht);

% Normal vector (inward)
n = zeros(3,nt,np);
for i = 1:nt
    for j = 1:np
        n(:,i,j) = cross(dxtt(:,i,j),dxpt(:,i,j));
        n(:,i,j) = -n(:,i,j)./norm(n(:,i,j));
    end
end
% Jacobian
J = sqrt(dot(dxtt,dxtt).*dot(dxpt,dxpt) - dot(dxtt,dxpt).^2);
J = squeeze(J);

end
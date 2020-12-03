function fnm = SpTNew(legs, f, wg, dphi,ord)
% New way to do spherical harmonic transforms. Assumes a uniform mesh...

% So we can handle a 3D transform
[nt,np,p] = size(legs);

if(np == 1)
    [nt,np] = size(f);
end

% In case we only want to go up to a certain order
if(nargin == 5)
    p = ord;
% Via quadratic formula    
else
    p = (-3+sqrt(9 - 4*(2-2*p)))/2; 
end

% For the FFTs in phi at constant thetas
gm = zeros(nt,np);

fnm = zeros((p+1)^2,1);

% Do the FFT, good for each constant value of theta
for i = 1:nt
    gm(i,:) = fft(f(i,:))*dphi;
end

% Now do the integrals individually across thetas
it = 0;
% Number of indices behind current order
ncnt = 0;
for n = 0:p
    ncnt = ncnt + n;
    for m = -n:n
        it = it+1;
%       Exploit symmetry        
        if(m>=0)
            legint = legs(:,1,ncnt + m + 1);
            gint = gm(:,1+m);
        else
            legint = legs(:,1,ncnt - m + 1)*(-1)^(m);
            gint = conj(gm(:,1-m));
        end
        fnm(it) = sum(wg.*legint.*gint);
    end
end
    
end

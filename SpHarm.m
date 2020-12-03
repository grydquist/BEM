function Ynm = SpHarm(n, tht, phi, facs)
%   Spherical harmonics of one order, all degrees
    [t,p] = size(tht);
    Pnm = legendre(n, cos(tht));
    Ynm = zeros(2*n+1, t, p);
    it = 0;
    for m = -n:n
        it = it+1;
        a = (2*n+1)/(4*pi);
        b = facs(n-abs(m)+1)/facs(n+abs(m)+1);
        C = sqrt(a*b);
        if n~=0
            Pm = squeeze(Pnm(abs(m)+1,:,:));
        else
            Pm = Pnm;
        end
        Ynm(it,:,:) = C*Pm.*exp(1i*abs(m)*phi);
        if m<0
            Ynm(it,:,:) = (-1)^abs(m)*conj(Ynm(it,:,:));
        end
    end
end
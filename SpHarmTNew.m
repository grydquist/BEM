function Ynm = SpHarmTNew(n, tht, phi)
% Spherical harmonics of all orders and degrees
Ynm = cell(n+1,1);

[t,p] = size(tht);
legs = myleg(n,cos(tht));

Ynm{1} = ones(1,t,p)*0.5*sqrt(1/pi);

for i = 1:n
    Ym = zeros(2*i+1, t, p);
    for m = i:-1:-i
        im = m + i + 1;
        if m<0
            Ym(im,:,:) = (-1)^(-m)*conj(Ym(im-2*m,:,:));
        else
            il = i*(i+1)/2 + m + 1;
            Ym(im,:,:) = legs(:,:,il).*exp(1i*m*phi);
        end
    end
    Ynm{i+1} = Ym;
end

end
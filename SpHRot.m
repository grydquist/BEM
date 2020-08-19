function f = SpHRot(fmn, a, b)

f = zeros(size(fmn));
p = sqrt(length(f))-1;

% Loop over harmonic order
it = 0;
for n = 0:p
    frng = it + 1:it + 2*n+1;
%   Loop over harmonic degree we're trying to calculate
    for mp = -n:n
        Dmm = zeros(1,2*n+1);
        im1 = 0;
        it = it+1;
%       Loop over harmonic degree we're using to calculate
        for m = -n:n
            Smm = 0;
            im1 = im1 + 1;
            for s = max(0,m - mp):min(n+m,n-mp)
                Smm = Smm + (-1)^s*(cos(b/2)^(2*(n-s)+m-mp)*sin(b/2)^(2*s-m+mp))/(factorial(n+m-s)*factorial(s)*factorial(mp-m+s)*factorial(n-mp-s));
            end 
            dmm = (-1)^(mp-m)*(factorial(n+mp)*factorial(n-mp)*factorial(n+m)*factorial(n-m))^0.5*Smm;
            Dmm(im1) = dmm*exp(1i*m*a);
        end
        f(it) = f(it) + sum(Dmm.*fmn(frng));
    end
end


end 
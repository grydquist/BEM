function f = SpHRot(fmn, a, b, c)

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
%             Dmm(im1) = exp(1i*m*a)*dmm*exp(1i*m*c);
            Dmm(im1) = exp(1i*m*a)*dmm;
        end
        f(it) = f(it) + sum(Dmm.*fmn(frng));
    end
end

% Here's the deal: the rotation was not working by doing it all at once. As
% it's written in the Sorgentone paper, the fact that we rotate around Z
% forward and backward makes the exp terms cancel. However, if you do
% things one at a time, it works. For example: start with f1, and rotate
% around Z (b = c = 0) to get f2. Rotate f2 around Y (a = c = 0) to get f3.
% Rotate f3 around Z again the negative amount (a = b = 0) to get the final
% f. It worksd doing the first Z and Y rotation together, but not the last
% one, so that is reflected here. There likely is a better way, but this
% works for now.

% Also, the rotations are in the right order, but JUST REMEMBER that the
% rotations are actually carried out on the vector in the order c->b->a


% Second loop over harmonic order to do last Z rotation.
it = 0;
for n = 0:p
%   Loop over harmonic degree we're trying to calculate
    for m = -n:n
        it = it+1;
        f(it) = f(it)*exp(1i*m*c);
    end    
end

end

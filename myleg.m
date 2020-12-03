function P = myleg(n, x)
% My own associated legendre function, already normalized. Computes values
% for ALL n & m (m >= 0) up to n.

% Size of grids we're evaluating at
[t, p] = size(x);

P=zeros(t,p,0.5*(n+1)*(n+2));

% Start with m,n == 0
P(:,:,1) = 1/sqrt(4*pi);
if(n == 0); return; end

% Get n = 1 values to start recurrence relationship
% m == 0 
P(:,:,2) = sqrt(3*0.25/pi)*x;
% m == 1
P(:,:,3) = -sqrt(3*0.125/pi)*(1 - x.^2).^(0.5);
if(n == 1); return; end

% Odd & even factorial to keep track of
facto = 1;
facte = 2;

it = 3;
for i = 2:n
    for m = 0:i
        it = it+1;
%       Analytical form for when m == n
        if(m == i)
            facto = facto*(2*m - 1);
            facte = facte*2*m*(2*m-1);
            P(:,:,it) = (-1)^m*sqrt((2*m + 1)/(4*pi*facte))*facto ...
                      * (1 - x.^2).^(0.5*m);
%       Recurrence relationship for m == n-1
        elseif(m == i-1)
            P(:,:,it) = sqrt(2*m + 3)*x.*P(:,:,it-i);
%      General recurrence relationship
        else
            P(:,:,it) = sqrt((4*i^2 - 1)/(i^2-m^2))*(x.*P(:,:,it - i) ...
                      - sqrt(((i - 1)^2 - m^2)/(4*(i-1)^2 - 1))*P(:,:,it - 2*i +1));
        end
    end
end
end
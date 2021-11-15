%% Periodic Stokeslet

function G = Gijp(xn, bxs, bv, eps1)
% xn: distance in original cell
% bxs: number of boxes to sum over
% bv: basis lattices vectors (columns of 3x3 matrix)
% eps1: parameter that weights importance of real and fourier parts
%       If not present use default below (usually balances 2 parts well)

G = zeros(3,3);

% Volume of boxes
tau = dot(cross(bv(:,1),bv(:,2)),bv(:,3));

% reciprocal basis vectors
kv = bv;
kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3));
kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1));
kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2));

if ~exist('eps1','var')
%   third parameter does not exist, so default it to something
    eps1 = pi^0.5/tau^(1/3);
end

%calculate volume of boxes

% Loop over all of the boxes we wish to sum over
for i = -bxs:bxs
    for j = -bxs:bxs
        for k = -bxs:bxs
%           Get current addition from box, add to separation
            rc = i*bv(:,1) + j*bv(:,2) + k*bv(:,3);
            xc = xn + rc;
%           Calculate contribution from real part
            G = G + HshortG(xc,eps1);
            
%           Current wave number in box
            kc = i*kv(:,1) + j*kv(:,2) + k*kv(:,3);
%           Fourier contribution (only if not in box)
            if(~(i==0 && j==0 && k==0))
                G = G + HspecG(kc, eps1)/tau*cos(dot(kc,xn));
            end

        end
    end
end

end



% Hasimotos
function G = HshortG(x,eps1)
r = norm(x);
er = r*eps1;
C = erfc(er) - 2/sqrt(pi)*er*exp(-er^2);
D = erfc(er) + 2/sqrt(pi)*er*exp(-er^2);
G = eye(3)*C/r + x*x'*D/r^3;
end

function G = HspecG(k,eps1)
kn = norm(k);
w = kn/eps1;
G = 8*pi/eps1^4*(1/w^4 + 1/4/w^2      )*(kn^2*eye(3) - k*k')*exp(-w^2/4);
end
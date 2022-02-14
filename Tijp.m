%% Periodic Stresslet
%%%%%%%%%% Doesn't seem to quite work in practice, used Marin instead
function T = Tijp(xn, bxs, bv, eps1, n)
% xn: distance in original cell
% xs: original x location for self interaction
% n: normal vec
% bxs: number of boxes to sum over
% bv: basis lattices vectors (columns of 3x3 matrix)
% eps1: parameter that weights importance of real and fourier parts
%       If not present use default below (usually balances 2 parts well)

T = zeros(3,3);

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
            T = T + BshortT(xc,eps1,n);
            
%           Current wave number in box
            kc = i*kv(:,1) + j*kv(:,2) + k*kv(:,3);
%           Fourier contribution (only if not in box)
            if(~(i==0 && j==0 && k==0))
                T = T + BspecT(kc, eps1,n)/tau*cos(dot(kc,xn));
            else
%               Self-interaction for zero net flow
%               Don't think needed as the box will move with dist. vel??
%                 T = T + selfT(xs,n)/tau;
            end

        end
    end
end

end

%%%% !!!!!!!!!!!!!!!!!!
%%%% T may need a self interaction term for the fourier part if k = 0 (this
%%%% will go to zero for the G part, but not necessarily for T).
%%%% I can't really tell though if this is true...
%%%% In fact I think we almost certainly don't need it. It isn't a
%%%% self-interaction term, it's different


% Beenakker's
function T = BshortT(x,eps1, n)
r = norm(x);
rh = x/r;
er = r*eps1;
xer = exp(-er^2);
C = -6*erfc(er)/r^2 - eps1/sqrt(pi)/r*(12 + 8*er^2 - 16*er^4)*xer;
D = 8*eps1^3/sqrt(pi)*(2-er^2)*xer; % Disagreement in lit if there's extra r here. Don't think there should be tho
rdotn = dot(rh,n);

T = C*(rh*rh'*rdotn) + D*(eye(3)*rdotn + n*rh' + rh*n');
end

function T = BspecT(k,eps1,n)
kn = norm(k);
w = kn/eps1;
kdotn = dot(k,n);
Q1 = -1i*(-2/kn^4*(k*k')*kdotn + 1/kn^2*(k*n' + n*k' + eye(3)*kdotn));
Q2 = 8 + 2*w^2 + w^4;
xer = exp(-w^2/4);

T = pi*Q1*Q2*xer;
end

function T = selfT(x,n)
T = 8*pi*dot(x,n);
end
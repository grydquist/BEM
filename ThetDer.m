function dY = ThetDer(Y, phi, n, m, ord)
% Recursively calculates the partial derivative of order ord w.r.t. theta
% Ym dhould be harmonics of all degree at order n in a cell array

Ym = Y{n+1};
dY = zeros(size(phi));
% Make sure phi is a mesh grid!

% If this is the bottom of the recursive function (i.e. first derivative)
if(ord == 1)
%       Two parts added together        
    if(m>-n)
        dY = dY ... 
            -0.5*sqrt((n+m)*(n-m+1))*exp( 1i*phi).*squeeze(Ym(m + n    , :, :));
    end
    if(m< n)
        dY =  dY ... 
            +0.5*sqrt((n-m)*(n+m+1)).*exp(-1i*phi).*squeeze(Ym(m + n + 2, :, :));
    end
else
% Else, recursive function
%   Two parts added together        
    if(m>-n)
        dY =  dY ... 
            -0.5*sqrt((n+m)*(n-m+1)).*exp( 1i*phi).*ThetDer(Y, phi, n, m-1, ord-1);
    end
    if(m< n)
        dY =  dY ... 
            +0.5*sqrt((n-m)*(n+m+1)).*exp(-1i*phi).*ThetDer(Y, phi, n, m+1, ord-1);
    end
end

end
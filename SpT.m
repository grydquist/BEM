function fnm = SpT(Ynm, f, tht, phi,ord)
% This should work, but we shouldn't have tht and phi be separate from the
% size of Ynm. It should be that phi = 2(p+1) and tht = p+1


% So we can handle a 3D transform
[sz,t,p] = size(f);

% If it really is just one thing to do at once
if(p == 1)
    p = t;
    t = sz;
    sz = 1;
end

if(nargin ==5)
    n = ord;
else
    n = length(Ynm)-1;
end
dphi = phi(1,2) - phi(1,1);

[xs,wg] = lgwt(t,-1,1);

fnm = zeros((n+1)^2,sz);
it = 0;
% Hacky but gets the job done
if(sz == 1)
    
    for i = 0:n
        Y = Ynm{i+1};
        for m = -i:i
            ind = m + i + 1;
            Yt = squeeze(Y(ind,:,:));
            it = it+1;
            for j = 1:t
                for k = 1:p
                    fnm(it) = fnm(it) + conj(Yt(j,k))*f(j,k)*wg(j)*dphi;
                end
            end
        end
    end
    
else
    fnm = zeros(sz,(n+1)^2);
    for i = 0:n
        Y = Ynm{i+1};
        for m = -i:i
            ind = m + i + 1;
            Yt = squeeze(Y(ind,:,:));
            it = it+1;
            for j = 1:t
                for k = 1:p
                    fnm(:,it) = fnm(:,it) + conj(Yt(j,k))*f(:,j,k)*wg(j)*dphi;
                end
            end
        end
    end
    
end

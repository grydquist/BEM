function x = SpHReconst(f,Y,p)
% Reconstructs field given spherical harmonic coefficients

% Max order to reconstruct to
if(nargin ==3)
    pp=p;
%   Useless now but could speed things up
    if(nargin==4)
        sz = 1;
    else
        sz = 1;
    end
else
    pp = length(Y)-1;
    sz = 1;
end
ih = 0;
it = 0;
for n = 0:pp
    ih = ih+1;
    Ycur = Y{ih};
    if (n==0)
        if sz == 1
            x = zeros(size(squeeze(Ycur)));
        else
            
        end
    end
    im = 0;
    for m = -n:n
        it = it+1;
        im = im+1;
        Ymcur = squeeze(Ycur(im,:,:));
        x = x + f(it)*Ymcur;
    end
end
end
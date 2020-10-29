function x = SpHReconst(f,Y,p)
% Reconstructs field given spherical harmonic coefficients
if(nargin ==3)
    pp=p;
else
    pp = length(Y)-1;
end
ih = 0;
it = 0;
for n = 0:pp
    ih = ih+1;
    Ycur = Y{ih};
    if (n==0); x = zeros(size(squeeze(Ycur)));end
    im = 0;
    for m = -n:n
        it = it+1;
        im = im+1;
        Ymcur = squeeze(Ycur(im,:,:));
        x = x + f(it)*Ymcur;
    end
end
end
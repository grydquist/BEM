function fnm = SpT(Ynm, f, tht, phi)
% This should work, but we shouldn't have tht and phi be separate from the
% size of Ynm. It should be that phi = 2(p+1) and tht = p+1

[t,p] = size(f);

n = length(Ynm)-1;
dphi = phi(1,2) - phi(1,1);

[xs,wg] = lgwt(t,-1,1);
% ws = zeros(1,length(wg));
% for i = 1:t
%     for j = 0:n
%         ws(i) = ws(i) + 2*sin(tht(i)/2)*legendreP(j,xs(i));
%     end
%     ws(i) = ws(i)*wg(i);
% end

fnm = zeros((n+1)^2,1);
it = 0;
ih = 0;
Yls = zeros((n+1)^2,(n+1)^2);
for i = 0:n
    Y = Ynm{i+1};
    for m = -i:i
        ih = ih+1;
        ind = m + i + 1;
        Yt = squeeze(Y(ind,:,:));
        it = it+1;
        ic = 0;
        for j = 1:t
            for k = 1:p
                ic = ic+1;
                fnm(it) = fnm(it) + conj(Yt(j,k))*f(j,k)*wg(j)*dphi;
                Yls(ih,ic) = Yt(j,k);
            end
        end
    end
end


function fnm = SpT(Ynm, f, tht, phi)
[t,p] = size(f);



n = length(Ynm)-1;
dphi = phi(2) - phi(1);

[xs,wg] = lgwt(t,-1,1);
ws = zeros(1,length(wg));
for i = 1:t
    for j = 0:n
        ws(i) = ws(i) + 2*sin(tht(i)/2)*legendreP(j,cos(tht(i)));
    end
    ws(i) = ws(i)*wg(i);
end

fnm = zeros((n+1)^2,1);
it = 0;
for i = 0:n
    Y = Ynm{i+1};
    for m = -i:i
        ind = m + i + 1;
        Yt = squeeze(Y(ind,:,:));
        it = it+1;
        for j = 1:t
            for k = 1:p
                fnm(it) = fnm(it) + Yt(j,k)*f(j,k)*ws(j)*dphi;
            end
        end
    end
end
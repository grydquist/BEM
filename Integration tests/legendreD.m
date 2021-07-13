function a = legendreD(k,x) 

if(k < 1)
    a = 0*x;
else
%   LegendreP is very slow b/c is uses symbolic
    b = legendre(k-1,x);
    c = legendre(k,x);
    a = k*(b(1,:) - x.*c(1,:))./(1-x.^2);
end

end
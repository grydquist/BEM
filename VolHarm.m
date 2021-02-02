function V = VolHarm(x,xt,xp,wg,tht,dphi)
%VOLHARM Calculates volume of spherical harmonics shape
V = 0;
for i = 1:length(tht)
    for j = 1:2*length(tht)
%       Integrand
        intgd = x(3,i,j)*(xt(1,i,j)*xp(2,i,j) - xp(1,i,j)*xt(2,i,j));
        
%       Integrate via GQ
        V = V + intgd*wg(i)*dphi/sin(tht(i));        
    end
end
end


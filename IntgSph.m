function intg = IntgSph(f, J, dphi, tht, wg)
[nt,np] = size(f);

intg = 0;
for i = 1:nt
    for j = 1:np
%         intg = intg + f(i,j)*wg(i)*J(i,j)*dphi/sin(tht(i));
        intg = intg + f(i,j)*wg(i)*J(i,j)*dphi;
    end
end

end
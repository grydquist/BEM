function SA = SAHarm(J,wg,tht,dphi)
%VOLHARM Calculates volume of spherical harmonics shape
SA = 0;
for i = 1:length(tht)
    for j = 1:2*length(tht)
%       Integrate via GQ
        SA = SA + J(i,j)*wg(i)*dphi/sin(tht(i));        
    end
end
end
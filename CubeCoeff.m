function [x1mn,x2mn,x3mn] = CubeCoeff(Y, tht, phi,ord)
nt = length(tht);
np = length(phi);
x1 = zeros(nt,np);
x2 = x1;
x3 = x1;

for i = 1:nt
    for j = 1:np
%         csphi = max([abs(cos(phi(j))), abs(sin(phi(j)))]);
%        if(tht(i) < pi/2 - atan2(1,csphi))
        if(tht(i)<pi/4)
            x1(i,j) = tan(tht(i))*cos(phi(j));%sqrt(2)*atan2(tht(i),1)*cos(phi(j));
            x2(i,j) = tan(tht(i))*sin(phi(j));%sqrt(2)*atan2(tht(i),1)*sin(phi(j));
            x3(i,j) = 1;
%         elseif(abs(tht(i) - pi) < pi/2 - atan2(1,csphi))
        elseif(tht(i) > 3*pi/4)
            x1(i,j) = tan(abs(tht(i)-pi))*cos(phi(j));%sqrt(2)*atan2(abs(tht(i) - pi),1)*cos(phi(j));
            x2(i,j) = tan(abs(tht(i)-pi))*sin(phi(j));%sqrt(2)*atan2(abs(tht(i) - pi),1)*sin(phi(j));
            x3(i,j) = -1;
        else
            if(tht(i) > pi/2)
                thtmp = pi/2 - abs(tht(i) - pi);
                x3(i,j) = -tan(thtmp);
            else
                thtmp = pi/2 - tht(i);
                x3(i,j) = tan(thtmp);
            end                
            x1(i,j) = cos(phi(j));
            x2(i,j) = sin(phi(j));
            x1(i,j) = sqrt(2)*cos(phi(j))*sin(tht(i));
            x2(i,j) = sqrt(2)*sin(phi(j))*sin(tht(i));
            x3(i,j) = sqrt(2)*cos(tht(i));
        end
    end
end

if(nargin ==4)
    n = ord;
else
    n = length(Y)-1;
end

x1mn = SpT(Y,x1,tht,phi,n);
x2mn = SpT(Y,x2,tht,phi,n);
x3mn = SpT(Y,x3,tht,phi,n);

end
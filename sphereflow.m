Ux = zeros(nx,ny);
Uy = zeros(nx,ny);

for i = 1:nx
    for j =1:ny
        theta = atan2(YY(j),XX(i));
        r = sqrt(XX(i)^2 + YY(j)^2);
        if r>=1
            ur = cos(theta)*(1 - 3/2/r + 1/2/r^3);
            ut = -sin(theta)*(1 - 3/4/r - 1/4/r^3);
            Ux(i,j) = ur*cos(theta) - ut*sin(theta);
            Uy(i,j) = ur*sin(theta) + ut*cos(theta);
        end
    end
end

quiver(XX,YY,Ux,Uy)
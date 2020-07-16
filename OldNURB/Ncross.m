function [c, cm, n] = Ncross(Xx, Xe)
% Evaluates the cross product of the derivatives at all xi's and eta's

ev1 = length(Xx(1,:,1));
ev2 = length(Xx(1,1,:));



c = zeros(3,ev1,ev2);
n = c;
cm =  zeros(ev1,ev2);

for i = 1:ev1
    for j = 1:ev2
        c(1,i,j) = Xx(2,i,j)*Xe(3,i,j) - Xx(3,i,j)*Xe(2,i,j);
        c(2,i,j) = Xx(3,i,j)*Xe(1,i,j) - Xx(1,i,j)*Xe(3,i,j);
        c(3,i,j) = Xx(1,i,j)*Xe(2,i,j) - Xx(2,i,j)*Xe(1,i,j);
        cm(i,j)  = sqrt(c(1,i,j)^2 + c(2,i,j)^2 + c(3,i,j)^2);
        n(1,i,j) = -c(1,i,j)/cm(i,j);
        n(2,i,j) = -c(2,i,j)/cm(i,j);
        n(3,i,j) = -c(3,i,j)/cm(i,j);
    end
end 

end
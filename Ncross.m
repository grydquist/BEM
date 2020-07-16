function [c, cm, n] = Ncross(Xx, Xe)
% Evaluates the cross product of the derivatives at all xi's and eta's

ev1 = length(Xx(1,:,1,1,1));
ev2 = length(Xx(1,1,:,1,1));
ev3 = length(Xx(1,1,1,:,1));
ev4 = length(Xx(1,1,1,1,:));

c = zeros(3,ev1,ev2,ev3,ev4);
n = c;
cm =  zeros(ev1,ev2,ev3,ev4);
for i = 1:ev1
    for j = 1:ev2
        for k = 1:ev3
            for l = 1:ev4
                c(1,i,j,k,l) = Xx(2,i,j,k,l)*Xe(3,i,j,k,l) - Xx(3,i,j,k,l)*Xe(2,i,j,k,l);
                c(2,i,j,k,l) = Xx(3,i,j,k,l)*Xe(1,i,j,k,l) - Xx(1,i,j,k,l)*Xe(3,i,j,k,l);
                c(3,i,j,k,l) = Xx(1,i,j,k,l)*Xe(2,i,j,k,l) - Xx(2,i,j,k,l)*Xe(1,i,j,k,l);
                cm(i,j,k,l)  = sqrt(c(1,i,j,k,l)^2 + c(2,i,j,k,l)^2 + c(3,i,j,k,l)^2);
                n(1,i,j,k,l) = -c(1,i,j,k,l)/cm(i,j,k,l);
                n(2,i,j,k,l) = -c(2,i,j,k,l)/cm(i,j,k,l);
                n(3,i,j,k,l) = -c(3,i,j,k,l)/cm(i,j,k,l);
            end
        end
    end
end 

end
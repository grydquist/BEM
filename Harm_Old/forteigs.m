function e = forteigs(A)
    p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3);

    I = eye(3);
    q = (A(1,1) + A(2,2) + A(3,3))/3;
    p2 = (A(1,1) - q)*(A(1,1) - q) + (A(2,2) - q)*(A(2,2) - q) ...
        + (A(3,3) - q)*(A(3,3) - q) + 2*p1;
    p = sqrt(p2/6);
    B = (1/p)*(A - q*I);
    r = det(B)*0.5;

    if(r < -1)
        phi = pi/3;
    elseif(r > 1) 
        phi = 0;
    else
        phi = acos(r)/3;
    end

    e(1) = q + 2*p*cos(phi);
    e(3) = q + 2*p*cos(phi + 2*pi/3);
    e(2) = 3*q - e(1) - e(3);
    
end
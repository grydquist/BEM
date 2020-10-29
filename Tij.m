function A = Tij(r,n)
    mri = 1/norm(r);
    A = diag(r.^2);
    A(1,2) = r(1)*r(2);
    A(1,3) = r(1)*r(3);
    A(2,3) = r(3)*r(2);
    A(2,1) = A(1,2);
    A(3,1) = A(1,3);
    A(3,2) = A(2,3);
    A = -6*A*mri^5*(r(1)*n(1) + r(2)*n(2) + r(3)*n(3));
end
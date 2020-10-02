function A = Gij(r)
    mri = 1/norm(r);
    A = diag(r.^2);
    A(1,2) = r(1)*r(2);
    A(1,3) = r(1)*r(3);
    A(2,3) = r(3)*r(2);
    A(2,1) = A(1,2);
    A(3,1) = A(1,3);
    A(3,2) = A(2,3);
    A = A*mri^3 + eye(3)*mri;
end
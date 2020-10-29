function Ynm = SpHarmT(n, tht, phi, myfacs)
% Spherical harmoonics of all orders and degrees
Ynm = cell(n+1,1);
for i = 0:n
    Ynm{i+1} = SpHarm(i,tht,phi,myfacs);
end
end
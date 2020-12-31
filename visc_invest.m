%https://www.sciencedirect.com/science/article/pii/S0096300319302292


phin = 10;
gamn = 100;
mus = zeros(gamn,phin);
gams = zeros(gamn,1);
sigs = mus;
phis = zeros(phin,1);
for i = 1:gamn
    for j = 1:phin
        phis(j) = 4*j/100;
        gams(i) = i;
        mus(i,j) = mu_calc(phis(j),gams(i));
        sigs(i,j) = mus(i,j)*gams(i);
    end
end

figure
hold on
for i = 1:phin
    plot(gams,mus(:,i)./mus(:,1))
end

figure
hold on
for i = 1:phin
    plot(sigs(:,i),mus(:,i)./mus(:,1))
end

figure
% loglog(0,0)
hold on
for i = 1:gamn
    plot(phis,mus(i,:))
end

function n0 = n0_calc(phi)
a1 = .1293;
a2 = -1.7912;
a3 = 7.6047;

n0 = 0.0736*(a1 + a2*phi + a3*phi^2);
end

function ninf = ninf_calc(phi)
b1 = .2114;
b2 = .9067;
b3 = 1.9879;

ninf = 0.005*(b1 + b2*phi + b3*phi^2);
end

function lam = lam_calc(phi)
c1 = 2.4024;
c2 = -2.381;

lam = 14.81*(c1 + c2*phi^2);
end

function mu = mu_calc(phi, gam)
ninf = ninf_calc(phi);
n0 = n0_calc(phi);
lam = lam_calc(phi);
mu = ninf + (n0 - ninf)*((1+log(1 + lam*gam))/(1 + lam*gam));
end


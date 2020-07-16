chi = linspace(-pi/2,pi/2,101);
alph = 1.38581894;
a = .72;
x = a*alph*sin(chi);
y = a*alph/2*(.207+2.003*sin(chi).^2-1.123*sin(chi).^4).*cos(chi);
plot(x,y)
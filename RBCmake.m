function [x,r] = RBCmake(tht, phi)
    nt = length(tht);
    np = length(phi);
    x = zeros(3,nt,np);
    r = zeros(nt,np);
    
    % First let's calculate half of the RBC for all tht at phi = 0, then
    % just rotate the points.
    
    % This honestly is not the best parameterization!!!!!!!!!!!
    % Works well enough but fix later
    b = 1;
    a = 1/1.1;
    ba = b/a;
    tt = tht - pi/2;
    
%     x0 = a*alph*sin(tht);
%     y0 = a*alph/2*(.207+2.003*sin(tht).^2-1.123*sin(tht).^4).*cos(tht);
%     r0 = sqrt(x0.^2+y0.^2);
    r0 = sqrt(a^2*cos(2*tt) + sqrt(ba^4 + sin(2*tt).^2));

    for i = 1:np
        x(1,:,i) = r0.*sin(tht)*cos(phi(i));
        x(2,:,i) = r0.*sin(tht)*sin(phi(i));
        x(3,:,i) = r0.*cos(tht);
        r(:,i) = r0;
    end

end
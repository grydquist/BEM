px = [-.5,-1.5];
S0 = 3;
pu = [0,0];
dt = 0.025;
St = 1;
maxt = 30;
mts = maxt/dt;
xs = zeros(mts,2);
ts = xs(:,1);
xa = [0,0];

for i = 1:mts
    S = rand*S0;
    vc = [px(1)*S,px(2)*-S];
    dv = 1/St*(vc - pu);
    pu = pu + dt*dv;
    px = px + dt*pu;
    
    if px(1) > 1
        px(1) = px(1) - 2;
        S0 = -S0;
        xa(1) = xa(1) + 2;
    end
    if px(2) > 1
        px(2) = px(2) - 2;
        S0 = -S0;
        xa(2) = xa(2) + 2;
    end
    if px(1) < -1
        px(1) = px(1) + 2;
        S0 = -S0;
        xa(1) = xa(1) - 2;
    end
    if px(2) < -1
        px(2) = px(2) + 2;
        S0 = -S0;
        xa(2) = xa(2) - 2;
    end
    xs(i,:) = px + xa;
    if i>1
        ts(i) = ts(i-1) + dt;
    end
end
plot(xs(:,1),xs(:,2))

%% 3d
px = [-.5,-1.5,-.5];
S0 = 3;
pu = [0,0,0];
dt = 0.025;
St = 1;
maxt = 30;
mts = maxt/dt;
xs = zeros(mts,3);
ts = xs(:,1);
xa = [0,0,0];
dU = zeros(3,3,mts);

for i = 1:mts
    S = rand*S0;
    vc = [px(1)*-0.5*S, px(2)*-0.5*S, px(3)*S];
    dv = 1/St*(vc - pu);
    pu = pu + dt*dv;
    px = px + dt*pu;
    
    if px(1) > 1
        px(1) = px(1) - 2;
        S0 = -S0;
        xa(1) = xa(1) + 2;
    end
    if px(2) > 1
        px(2) = px(2) - 2;
        S0 = -S0;
        xa(2) = xa(2) + 2;
    end
    if px(3) > 1
        px(3) = px(3) - 2;
        S0 = -S0;
        xa(3) = xa(3) + 2;
    end
    
    if px(1) < -1
        px(1) = px(1) + 2;
        S0 = -S0;
        xa(1) = xa(1) - 2;
    end
    if px(2) < -1
        px(2) = px(2) + 2;
        S0 = -S0;
        xa(2) = xa(2) - 2;
    end
    if px(3) < -1
        px(3) = px(3) + 2;
        S0 = -S0;
        xa(3) = xa(3) - 2;
    end
    xs(i,:) = px + xa;
    if i>1
        ts(i) = ts(i-1) + dt;
    end
    
    dU(1,1,i) = -0.5*S;
    dU(2,2,i) = -0.5*S;
    dU(3,3,i) = S;
end
plot3(xs(:,1),xs(:,2),xs(:,3))

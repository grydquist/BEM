% !!!!!!! CHECK TO SEE RBC RETURNS TO OG SHAPE
%% Given a deformed and reference configuration, calculate the traction on
% the surface of the object, for an infinitesimal 2d membrane of spherical
% harmonic topology.

% Start with some matrial constants.
B = 0.005;
C = 100;
Eb = 1;%0.0669; % Non-dimensional

%% Defining the shapes

% Order of harmonics
p = 3;

% Our grid to evaluate integrals over the surface where we need these forces.
gf = 1;
np = gf*2*(p+1);
dphi = 2*pi/np;
phi = 0:dphi:dphi*(np-1)';
nt = gf*(p+1);
[xs,wg] = lgwt(nt,-1,1);
tht = acos(xs);
[ph,th] = meshgrid(phi,tht);

% Reference
xmnR = 2*sqrt(pi);

% Order of geometric parameterization
pxmn = sqrt(length(xmnR))-1;

% Deformed
xmn = 2*sqrt(pi);

% Harmonics evaluated at the given thetas/phis
Yt = SpHarmT(pxmn,th,ph);

% Derivatives and whatnot in reference and deformed
[xR, nR, JR, xdtR, xdpR, kR, kdR] = SpHDer(xmnR, Yt, th,ph);
[x , n , J , xdt , xdp , k , kd ] = SpHDer(xmn , Yt, th,ph);
tauab = zeros(2,2);
tau11 = zeros(nt,np);
tau12 = tau11;
tau21 = tau11;
tau22 = tau11;
q1 = tau11;
q2 = tau11;



%       Bending moment !!!!!!!!!!! NOT THE MOST ACCURATE WAY OF DOING IT
%       And honestly probably not the best either as far as derivatives go.
%       Instead maybe I should try Zhao's
        q = Eb*SpBend(xmn,th,ph, Yt, xdt,xdp,kR,kdR);

% Left Cauchy-Green tensor/principal stretches/dirs
for i = 1:nt
    for j = 1:np
        
%       Contravariant form of reference tangent
        cdtR = cross(xdpR(:,i,j),-nR(:,i,j))/JR(i,j);
        cdpR = cross(-nR(:,i,j),xdtR(:,i,j))/JR(i,j);
        
%       Deformation gradient tensor
        F = xdt(:,i,j)*cdtR' + xdp(:,i,j)*cdpR';
        
%       Surface projection operator
        P = eye(3) - n(:,i,j)*n(:,i,j)';
        
%       Surface deformation gradient tensor
        Fs = P*F;
        
%       Cauchy-Green
        V2 = F*F';
        
%       Principal strains
        [ev,lams] = eigs(V2);
        es = sqrt(diag(lams));
        
%       Strain invariants
        I1 = es(1)^2 + es(2)^2 - 2;
        I2 = es(1)^2*es(2)^2 - 1;
        
%       In plane tension
        tau = B/(2*es(1)*es(2))*(I1 + 1)*V2 ...
            + 0.5*es(1)*es(2)*(C*I2-B)*P;
        
%       Principal tensions
        taup1 = es(1)/es(2)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);
        taup2 = es(2)/es(1)*(B/2*(2*I1+1)) + es(1)*es(2)*(-B/2+C/2*I2);

%       Contra- and covariant bases for deformed
        a1 = xdt(:,i,j);
        a2 = xdp(:,i,j);
        c1 = cross(a2,-n(:,i,j))/J(i,j);
        c2 = cross(-n(:,i,j),a1)/J(i,j);
        
%       Matrix in surface coordinates (contravariant)
        tauab(1,1) = c1'*tau*c1;
        tauab(1,2) = c1'*tau*c2;
        tauab(2,1) = c2'*tau*c1;
        tauab(2,2) = c2'*tau*c2;
        
%       Put them in a form amenable to a spherical harmonic transform
        tau11(i,j) = tauab(1,1);
        tau12(i,j) = tauab(1,2);
        tau21(i,j) = tauab(2,1);
        tau22(i,j) = tauab(2,2);
        q1(i,j) = q(1,i,j);
        q2(i,j) = q(2,i,j);
        
    end
end

Yp = SpHarmT(p,th,ph);
f11 = SpT(Yp,tau11,th,ph);
f12 = SpT(Yp,tau12,th,ph);
f21 = SpT(Yp,tau21,th,ph);
f22 = SpT(Yp,tau22,th,ph);
fq1 = SpT(Yp,q1,th,ph);
fq2 = SpT(Yp,q2,th,ph);

% Package these up to send into ForceCalc
ftot = [f11,f12,f21,f22,fq1,fq2];


% Now we can calculate the forces by taking the curvature tensors and
% covariant derivatives.
f = real(ForceCalc(xmn,ftot,Yp,th,ph, xdt, xdp));


%% What's wrong with this program:
% The functions are all over the place and I'm doing the same thing over
    % and over again, like getting the Christoffel symbols and whatnot
    % (e.g. SpBend and ForceCalc)
% So maybe instead I can make these not be functions, but instead be
    % scripts so everyone can have access to the variables.


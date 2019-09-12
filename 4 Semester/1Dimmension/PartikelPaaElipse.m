clear all; close all;clc;
numberOfPoints = 400;
NumberOfPlots = 2

syms t u v
% r = [cos(t);sin(t);z]
parameter = linspace(1.1,10,5)
for ind = 1:length(parameter)
x = 100*cos(u);
y = parameter(ind)*100*sin(u);
z = v;

r(u,v) = [x;y;z];
nor =@(x) sqrt(x'*x);

ru = simplify(diff(r,u));
ruu = simplify(diff(ru,u));
rv = simplify(diff(r,v));
rvv = simplify(diff(rv,v));
ruv = simplify(diff(rv,u));
n = simplify(cross(rv,ru)/nor(cross(rv,ru)));

%% First and second fundemental form
L = simplify(ruu'*n);
M = simplify(ruv'*n);
N = simplify(rvv'*n);
E = simplify(-ru'*ru);
F = simplify(-rv'*ru);
G = simplify(-rv'*rv);

%% Principal radii of curvature
H = simplify((E*N+G*L-2*F*M)/(2*E*G-2*F^2));
K = simplify((L*N-M^2)/(E*G-F^2));

k1 = simplify(H + sqrt(H^2-K));
k2 = simplify(H - sqrt(H^2-K));

%% Effective potential
U = simplify(1/4*(k1-k2)^2);
% Ufun = matlabFunction(U+eps^3*(u+v))
%% Simulation and plot
% N = NumberOfPoints; %number of simulated points minus 2 

g = [-E,-F;-F,-G];
us = linspace(-pi,pi,numberOfPoints+1); % Her tages det sidste punkt med til plots. Det gør us n+1 lang

Ufun = matlabFunction(U);
Us = Ufun(us(:)).*ones(size(us(:)));

% figure
rFun = matlabFunction(r);
rs = rFun(us,zeros(size(us)));
xs = rs(1,:);
ys = rs(2,:);
% plot(xs,ys)
% axis equal
% title('The elipse') % Titel på plot
% xlabel('x') % Navn på x-akse
% ylabel('y') % Navn på y-akse
% hold on 
% plot3(xs,ys,Us)

% figure
% title('"Energy"') % Titel på plot
% xlabel('theta') % Navn på x-akse
% ylabel('U') % Navn på y-akse
% plot(us,Us)

n = numberOfPoints;
M = sparse(zeros(n,n));

du = us(2)-us(1);

g1 = det(g);

% g = matlabFunction(g+eps^3*(u));
gdu = matlabFunction(simplify(diff(g1,u))+eps^3*(u));
gI = matlabFunction(g1^-1+eps^3*(u));
% dgIu = matlabFunction(simplify(diff(gI,u))+eps^3*(u));
for K = 2:(n-1)
    gNum = gI(us(K));
    guNum = gdu(us(K));

    
    M(K,K-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M(K,K)   = 2*gNum/(du^2)-Us(K);
    M(K,K+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;
end

gNum = gI(us(1));
guNum = gdu(us(1));
M(1,n)   =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(1,1)   = 2*gNum/(du^2)-Us(1);
M(1,1+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;
    
gNum = gI(us(n));
guNum = gdu(us(n));
M(n,n-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(n,n)   = 2*gNum/(du^2)-Us(n);
M(n,1)   =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;

% figure
% surf(M)
opts.tol = 1e-14
[V,D] = eigs(M,NumberOfPlots,'sr',opts);
[E,I] = sort(diag(D));

%% plot
for i =1:1
    figure
    psi2 = conj(V(:,I(i))).*V(:,I(i));
    psii = [psi2;psi2(1)];

    plot(xs,ys)
    hold on
    plot3(xs,ys,psii)

    title(['Energy nr: ' num2str(i) ' of E= ' num2str(E(i))])
    xlabel('x')
    ylabel('y')
    figure
    plot(us,psii)
end
Esave(ind,:) = E'
end

% M = zeros(n-1);
% M(1,n-1) = -gf(theta(n-1))/dT^2;
% M(1,1) = gf(theta(1))*2/dT^2+3/2*dgf(theta(1))/dT-uf(theta(1));
% M(1,2) = -gf(theta(2))/dT^2-3/2*dgf(theta(2))/dT;
% 
% M(n-1,n-2) = -gf(theta(n-2))/dT^2;
% M(n-1,n-1) = gf(theta(n-1))*2/dT^2+3/2*dgf(theta(n-1))/dT-uf(theta(n-1));
% M(n-1,1) = -gf(theta(1))/dT^2-3/2*dgf(theta(1))/dT;
% 
% 
% M;
% [V,D] = eig(M);
% [E,I] = sort(diag(D));
% E(1:10)
% I(1:10)
% % [E,I] = min(min(D))
% 
% psi1 = [V(:,I(2));V(1,I(2))];
% s = cumtrapz(sqrt(gf(theta)).*dT);
% K = sqrt(trapz(s,psi1.^2))
% psi1=psi1/K;
% figure
% plot(s,psi1)
% title(['First wave function, energy' num2str(E(2))]) % Titel på plot
% xlabel('s') % Navn på x-akse
% ylabel('psi') % Navn på y-akse
% 
% psi2 = [V(:,I(1));V(1,I(1))];
% K = sqrt(trapz(s,psi2.^2))
% 
% psi2=psi2./K;
% figure
% plot(s,psi2)
% title(['Second wave function, energy' num2str(E(1))]) % Titel på plot
% xlabel('s') % Navn på x-akse
% ylabel('psi') % Navn på y-akse
% 
% 
% %% Tilstand på en side
% figure
% title(['Combiend wave function, energy' num2str((E(1)+E(2))/2)]) % Titel på plot
% xlabel('s') % Navn på x-akse
% ylabel('psi') % Navn på y-akse
% plot(s,(psi2-psi1)/sqrt(2))
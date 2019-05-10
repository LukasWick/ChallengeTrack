clear all; close all;clc;
numberOfPoints = 500;
NumberOfPlots = 10;

syms t u v
% r = [cos(t);sin(t);z]
x = u;
y = u*eps^3;
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
us = linspace(-pi,pi,numberOfPoints+2); % Ekstra punkt i hver ende
Ufun = matlabFunction(U);
Us = Ufun(us(:)).*ones(size(us(:)));

figure
rFun = matlabFunction(r);
rs = rFun(us,zeros(size(us)));
xs = rs(1,:);
ys = rs(2,:);
plot(xs,ys)
axis equal
title('The elipse') % Titel p� plot
xlabel('x') % Navn p� x-akse
ylabel('y') % Navn p� y-akse
hold on 
plot3(xs,ys,Us)

figure
title('"Energy"') % Titel p� plot
xlabel('theta') % Navn p� x-akse
ylabel('U') % Navn p� y-akse
plot(us,Us)

n = numberOfPoints;
M = sparse(zeros(n,n));

du = us(2)-us(1);

g1 = det(g);

% g = matlabFunction(g+eps^3*(u));
gdu = matlabFunction(simplify(diff(g1,u))+eps^3*(u));
gI = matlabFunction(g1^-1+eps^3*(u));
% dgIu = matlabFunction(simplify(diff(gI,u))+eps^3*(u));
for K = 2:(n-1)
    gNum = gI(us(K+1));
    guNum = gdu(us(K+1));

    
    M(K,K-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M(K,K)   = 2*gNum/(du^2)-Us(K+1);
    M(K,K+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;
end

gNum = gI(us(1+1));
guNum = gdu(us(1+1));
% M(1,n)   =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(1,1)   = 2*gNum/(du^2)-Us(1+1);
M(1,1+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;
    
gNum = gI(us(n+1));
guNum = gdu(us(n+1));
M(n,n-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(n,n)   = 2*gNum/(du^2)-Us(n+1);
% M(n,1)   =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;

% figure
% surf(M)
[V,D] = eigs(M,NumberOfPlots,'sr');
[E,I] = sort(diag(D));

%% plot
for i =1:3
    figure
    psi2 = conj(V(:,I(i))).*V(:,I(i));
    psii = [0;psi2;0];

    plot(xs,ys)
    hold on
    plot3(xs,ys,psii)

    title(['Energy nr: ' num2str(i) ' of E= ' num2str(E(i))])
    xlabel('x')
    ylabel('y')
    figure
    plot(us,psii)
end



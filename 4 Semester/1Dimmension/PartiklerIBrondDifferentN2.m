clear all; close all;clc;

syms t u v
% r = [cos(t);sin(t);z]
x = u;
y = u;
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

Ufun = matlabFunction(U);

g1 = det(g);

gdu = matlabFunction(simplify(diff(g1,u))+eps^3*(u));
gI = matlabFunction(g1^-1+eps^3*(u));
rFun = matlabFunction(r);

numberOfPointsS = linspace(100,1000,40)


for index = 1:length(numberOfPointsS)
numberOfPoints = round(numberOfPointsS(index));
NumberOfPlots = 2

width = pi/sqrt(2)
us = linspace(-width,width,numberOfPoints+2); % Ekstra punkt i hver ende
Us = Ufun(us).*ones(size(us));

rs = rFun(us,zeros(size(us)));
xs = rs(1,:);
ys = rs(2,:);

n = numberOfPoints;
M = sparse(zeros(n,n));

du = us(2)-us(1);


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

opts.maxit = 1e+16

[V,D] = eigs(M,NumberOfPlots,'sr',opts);
[E,I] = sort(diag(D));

Energies(index) = E(1)
end
%%
figure
plot(numberOfPointsS,(Energies-1/4)*10^5 )
hold on
title('Error convergens infinite well')
set(gca,'FontSize',14) 
xlabel('Number of points')
ylabel('Error: \DeltaE  [10^{- 5}]')
xlim([0,1000])
ylim([-2.25,0.05])
grid on

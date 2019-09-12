clear all; close all;clc;
numberOfPoints =50;
NumberOfPlots = 3;

syms t u v
x = cos(u);
y = sin(u);
z = v;

% x = cos(u)*(0.2*cos(4*u)+2);
% y = sin(u)*(0.2*cos(4*u)+2);
% z = v;


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

figure
rFun = matlabFunction(r);
rs = rFun(us,zeros(size(us)));
xs = rs(1,:);
ys = rs(2,:);


%% Plotting the setup
plot(xs,ys)
axis equal
title('The elipse') % Titel på plot
xlabel('x') % Navn på x-akse
ylabel('y') % Navn på y-akse
hold on 
plot3(xs,ys,Us)

figure
title('"Energy"') % Titel på plot
xlabel('theta') % Navn på x-akse
ylabel('U') % Navn på y-akse
plot(us,Us)

n = numberOfPoints;

du = us(2)-us(1);

gdet = det(g);

dgdu = matlabFunction(simplify(diff(gdet,u))+eps^3*(u));
gInv = matlabFunction(gdet^-1+eps^3*(u));
M = sparse(zeros(n,n));

gInvNum = gInv(us);
dgduNum = dgdu(us);


for K = 2:(n-1)
    gNum  = gInvNum(K);
    guNum = dgduNum(K);
    M(K,K-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M(K,K)   = 2*gNum/(du^2);
    M(K,K+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;
end

gNum  = gInvNum(1);
guNum = dgduNum(1);
M(1,n)   =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(1,1)   = 2*gNum/(du^2);
M(1,1+1) =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;

gNum  = gInvNum(n);
guNum = dgduNum(n);
M(n,n-1) =  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum; 
M(n,n)   = 2*gNum/(du^2);
M(n,1)   =  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;

M2 = sparse(zeros(n^2,n^2));

interactionFun = @(L)  10./(L+0.01);

% interactionFun = @(L)  -exp(-L.^2);
% interactionFun = @(L) 0
%%
for i = (0:n:n^2-n)
    Ind1 = 1:n;
    Ind2 = i/n+1;
    L = sqrt((xs(Ind1)-xs(Ind2)).^2+(ys(Ind1)-ys(Ind2)).^2);
    Length(Ind2,:) = L;
    vMat = diag(interactionFun(L));
    Umat = diag(ones(n,1).*(Us(Ind1)+Us(Ind2)));
    M2(i+1:(i+n),i+1:(i+n)) = M - Umat + vMat;
end

%%

for L = 2:(n-1)
   gNum  = gInvNum(L);
   guNum = dgduNum(L);
   for K = 1:n 
    index = K + (L-1)*n;
    M2(index,index-n) =M2(index,index-n)  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M2(index,index)   =M2(index,index)+ 2*gNum/(du^2);
    M2(index,index+n) =M2(index,index+n)  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;  
   end
end
gNum  = gInvNum(1);
guNum = dgduNum(1);
for K = 1:n 
    index = K + (1-1)*n;
    M2(index,index-n+n^2) =M2(index,index-n+n^2)  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M2(index,index)   =M2(index,index)+ 2*gNum/(du^2);
    M2(index,index+n) =M2(index,index+n)  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;  
end
gNum  = gInvNum(n);
guNum = dgduNum(n);
for K = 1:n 
    index = K + (n-1)*n;
    M2(index,index-n) =M2(index,index-n)  -gNum/(du^2) - 1/(4*du)*gNum.^2*guNum;
    M2(index,index)   =M2(index,index)+ 2*gNum/(du^2);
    M2(index,index+n-n^2) =M2(index,index+n-n^2)  -gNum/(du^2) + 1/(4*du)*gNum.^2*guNum;  
end
%%

[V,D] = eigs(M2,NumberOfPlots,'sr');
[E,I] = sort(diag(D));

%%
% s = cumtrapz(sqrt((gInvNum(1:n)).^(-1)).*du);
for i = 1:NumberOfPlots
% psiNorm = conj(V(:,I(i))).*V(:,I(i));
% psi = reshape(psiNorm,[n,n]);
psi = reshape(V(:,I(i)),[n,n]);
psi1 = psi.^2*(sqrt((gInvNum(1:n)).^(-1)).*du)';
% psi1 = sum(psi')
psi1 = [psi1;psi1(1)];

psi2 = psi.^2'*(sqrt((gInvNum(1:n)).^(-1)).*du)';
psi2 = [psi2;psi2(1)];

figure
surf(psi)

figure
plot(xs,ys)
title('Solution') % Titel på plot
xlabel('x') % Navn på x-akse
ylabel('y') % Navn på y-akse
hold on 
plot3(xs,ys,psi1,'r','linewidth',1.5)
plot3(xs,ys,psi2,'b--','linewidth',1.5)
daspect([1,1,0.01])
end

%% Liniar kombination
psi = reshape(sqrt(2).*(V(:,I(1))+V(:,I(2))),[n,n]);
psi1 = psi.^2*(sqrt((gInvNum(1:n)).^(-1)).*du)';
% psi1 = sum(psi')
psi1 = [psi1;psi1(1)];

psi2 = psi.^2'*(sqrt((gInvNum(1:n)).^(-1)).*du)';
psi2 = [psi2;psi2(1)];

figure
surf(psi)

figure
plot(xs,ys)
title('Solution') % Titel på plot
xlabel('x') % Navn på x-akse
ylabel('y') % Navn på y-akse
hold on 
plot3(xs,ys,psi1,'r','linewidth',1.5)
plot3(xs,ys,psi2,'b--','linewidth',1.5)
daspect([1,1,0.01])

%% The "potential"

figure
[Us1,Us2] = meshgrid(us(1:end-1),us(1:end-1));
UU1 = -Ufun(Us1);
UU2 = -Ufun(Us2);

figure
surf(Us1,Us2,UU2+UU1)
figure
surf(Us1,Us2,interactionFun(Length))
figure
surf(Us1,Us2,UU2+UU1+interactionFun(Length))

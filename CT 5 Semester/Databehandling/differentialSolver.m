% clc; clear all; close all;
alpha = 0

T = 0.2500;
n = 60;
times = linspace(0,T,n);
dat0 = [0;2.3];

acc_n = 5
acc_points = normrnd(0,3,[acc_n,1]);
acc_points_x =linspace(0,T*dat0(2),acc_n)%[-0.2;-0.1;rand([3,1])*T*dat0(2)/2+T*dat0(2)/2,;T*dat0(2)+[0.1;0.2]];
accFun = @(x) spline(acc_points_x,acc_points,x);


[t,dat] = ode45(@(t,dat) diffFun(t,dat,accFun,alpha),times,dat0);
xs = dat(:,1);
vs = dat(:,2);

dat0 = [xs(end);-vs(end)];

[temp,dat] = ode45(@(t,dat) diffFun(t,dat,accFun,alpha),times,dat0);
t = [t;t(end)+temp]
xs = [xs;dat(:,1)];
vs = [vs;dat(:,2)];

%%
figure
plot(t,xs)
figure 
plot(t,vs)
figure
plot(xs,accFun(xs))
hold on 
plot(acc_points_x,acc_points,'.')
%%

function dDat = diffFun(t,dat,accFun,alpha)
x = dat(1);
v = dat(2);
dxdt = v;
dvdt = accFun(x)-alpha*v;
dDat = [dxdt;dvdt];
end


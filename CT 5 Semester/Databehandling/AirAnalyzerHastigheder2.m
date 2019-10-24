clc; clear all; close all;

n_poly=5;

P=@(a,x)a(1);
dP=@(a,x) 0;
ddP=@(a,x) 0;
for i=1:n_poly
P=@(a,x) P(a,x)+a(i+1)*x.^i;
dP=@(a,x) dP(a,x)+i*a(i+1)*x.^(i-1);
end
for i=2:n_poly
ddP=@(a,x) ddP(a,x)+i*(i-1)*a(i+1)*x.^(i-2);
end

% A=load('driftmåling1.txt');
A=load('drag.txt');

t1=A(:,1);
x=A(:,2);
y=A(:,3);

t=t1/((t1(end)-t1(1))/(length(t1)-1)*240);

% plot(x(1e4:end),'.')
% hold on
% plot(y(1e4:end),'.')

plot(x,'.')
hold on
plot(y,'.')


%%
% close all;
i1=[1 63];
i2=[60 140];

am=[];
X=[];
Y=[];
Ax=[];
Ay=[];
for i=1:2
ts=t(i1(i):i2(i));

xs=x(i1(i):i2(i));
ys=y(i1(i):i2(i));
figure
plot(xs,'.')
hold on
plot(ys,'.')

% b1=nlinfit(ts,xs,@(b,x)b(1)+b(2)*x,[1 1]);
% xs=xs-b1(1)-b1(2)*ts;
% b2=nlinfit(ts,ys,@(b,x)b(1)+b(2)*x,[1 1]);
% ys=ys-b2(1)-b2(2)*ts;

vx=gradient(xs,ts);
vy=gradient(ys,ts);
ax=gradient(vx,ts);
ay=gradient(vy,ts);

t0=ts(1);
a1(i,:)=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
a2(i,:)=nlinfit(ts,ys,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
figure
plot(ts,xs,'.')
hold on
plot(ts,P(a1(i,:),ts-t0))
figure
plot(ts,ys,'.')
hold on
plot(ts,P(a2(i,:),ts-t0))
end

t1=t(i1(1):i2(1))-t(i1(1));
for i = (i1(1):i2(1))-i1(1)+1
   t2(i) = fzero(@(t2) P(a1(1,:),t1(i))-P(a1(2,:),t2),t1(max(1,length(t1)-i)));     
end
t2 = t2';
Alpha = -(ddP(a1(2,:),t2)-ddP(a1(1,:),t1))./(dP(a1(2,:),t2)-dP(a1(1,:),t1));
figure
title('Alpha')
plot(Alpha,'.')
figure
plot(vx,Alpha,'.')


clc; clear all; close all;

n_poly=7;

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
% 
% A=load('driftmåling1.txt');
% t=A(:,1);
% x=A(:,2);
% y=A(:,3);
% 
% A=load('driftmåling2.txt');
% t=[t;A(:,1)];
% x=[x;A(:,2)];
% y=[y;A(:,3)];

A=load('drag.txt');
t=A(:,1);
x=A(:,2);
y=A(:,3);


plot(1:length(x),x)
figure 
% plot(x(1e4:end),'.')
% hold on
% plot(y(1e4:end),'.')
plot(x,'.')
hold on
plot(y,'.')


%%
% close all;

i1 = [62]
i2 = [140]

am=[];
X=[];
Y=[];
Ax=[];
Ay=[];
% for i=1:length(i1)
% for i=length(i1)-4:length(i1)
for i=1:length(i1)
ts=t(i1(i):i2(i));

xs=x(i1(i):i2(i));
ys=y(i1(i):i2(i));
% plot(xs,'.')
% hold on
% plot(ys,'.')

% b1=nlinfit(ts,xs,@(b,x)b(1)+b(2)*x,[1 1]);
% xs=xs-b1(1)-b1(2)*ts;
% b2=nlinfit(ts,ys,@(b,x)b(1)+b(2)*x,[1 1]);
% ys=ys-b2(1)-b2(2)*ts;

vx=gradient(xs,ts);
vy=gradient(ys,ts);
ax=gradient(vx,ts);
ay=gradient(vy,ts);

t0=ts(1);
a1=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1),statset('MaxIter',5000,'TolFun',1e-15,'TolX',1e-15));
a1=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
a2=nlinfit(ts,ys,@(a,t)P(a,t-t0),zeros(1,n_poly+1));

figure
plot(ts,xs,'.')
hold on
plot(ts,P(a1,ts-t0))
title(['x (',num2str(i),')'])


% 
% figure
% plot(ts,vx,'.')
% hold on
% plot(ts,dP(a1,ts-t0))
% title(['vx (',num2str(i),')'])
% 
% figure
% plot(ts,ax,'.')
% hold on
% plot(ts,ddP(a1,ts-t0))
% title(['ax (',num2str(i),')'])

% figure
% plot(ts,ys,'.')
% hold on
% plot(ts,P(a2,ts-t0))
% title(['y (',num2str(i),')'])
% 
% figure
% plot(ts,vy,'.')
% hold on
% plot(ts,dP(a2,ts-t0))
% title(['vy (',num2str(i),')'])
% 
% figure
% plot(ts,ay,'.')
% hold on
% plot(ts,ddP(a2,ts-t0))
% title(['ay (',num2str(i),')'])

time=linspace(t0,ts(end),1000);
am=[am max(sqrt(ddP(a1,time-t0).^2+ddP(a2,time-t0).^2))];

time=linspace(t0,ts(end),20);
X=[X P(a1,time-t0)];
Y=[Y P(a2,time-t0)];
Ax=[Ax ddP(a1,time-t0)];
Ay=[Ay ddP(a2,time-t0)];
end
am=max(am)

%%
% close all
figure
scale=min(sqrt((X(2:end)-X(1:end-1)).^2+((Y(2:end)-Y(1:end-1)).^2)))
ColorFieldV(X,Y,Ax,Ay,scale)
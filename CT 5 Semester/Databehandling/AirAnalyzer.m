clc; clear all; close all;

n_poly=7;
K=5;

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

A=load('...\Data\Drift1.txt');
t1=A(:,1);
x=A(:,2);
y=A(:,3);

t1=t1/((t1(end)-t1(1))/(length(t1)-1)*240);

A=load('...\Data\Drift2.txt');
t2=[A(:,1)];
x=[x;A(:,2)];
y=[y;A(:,3)];

t2=t2/((t2(end)-t2(1))/(length(t2)-1)*240);

t=[t1;t2];
dt = 1/240;
t = dt.*(1:length(t));
t = t';
% plot(x(1e4:end),'.')
% hold on
% plot(y(1e4:end),'.')
plot(x,'.')
hold on
plot(y,'.')


%%
close all;
% i1=[3 820 1605 2395 3108 3871 4552 4782 5698 6205 7153 8801 9542 10865 11864 12076 12982 13187];
% i2=[210 957 1744 2539 3257 4005 4650 4923 5817 6324 7277 8897 9637 10972 11954 12199 13070 13320];
% i1=[3 820 1605 2395 3108 3871 4552 4782 5698 6205 7153 8801 9542 10865 11864 12076 13187];
% i2=[210 957 1744 2539 3257 4005 4650 4923 5817 6324 7277 8897 9637 10972 11954 12199 13320];

i1=[3 820 1605 2395 3108 3871 4552 4782 5698 6205 7153 8801 9542 10865 12076 13187];
i2=[210 957 1744 2539 3257 4005 4650 4923 5817 6324 7277 8897 9637 10972 12199 13320];

am=[];
X=[];
Y=[];
Ax=[];
Ay=[];
Xs=[];
Ys=[];
for i=1:length(i1)-4
% for i=length(i1)-4:length(i1)
%     for i=1:1
ts=t(i1(i):i2(i));

xs=x(i1(i):i2(i));
ys=y(i1(i):i2(i));
Xs=[Xs;xs];
Ys=[Ys;ys];
% plot(xs,ys,'.')
% hold on
% plot(ys,'.')
% 
% b1=nlinfit(ts,xs,@(b,x)b(1)+b(2)*x,[1 1]);
% xs=xs-b1(1)-b1(2)*ts;
% b2=nlinfit(ts,ys,@(b,x)b(1)+b(2)*x,[1 1]);
% ys=ys-b2(1)-b2(2)*ts;

vx=gradient(xs,ts);
vy=gradient(ys,ts);
ax=gradient(vx,ts);
ay=gradient(vy,ts);

% dt = ts(3)-ts(4)
% vx=gradient(xs,dt);
% vy=gradient(ys,dt);
% ax=gradient(vx,dt);
% ay=gradient(vy,dt);
% figure
% plot(diff(ts))
% figure

t0=ts(1);
% a1=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1),statset('MaxIter',5000,'TolFun',1e-15,'TolX',1e-15));
a1=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
a2=nlinfit(ts,ys,@(a,t)P(a,t-t0),zeros(1,n_poly+1));

% figure 
% histogram(xs-P(a1,ts-t0),25)

% plotmyfigures(ts,xs,vx,ax,ys,vy,ay,P,dP,ddP,a1,a2,t0,i)

time=linspace(t0,ts(end),1000);
am=[am sqrt(ddP(a1,time(K:end-K)-t0).^2+ddP(a2,time(K:end-K)-t0).^2)];

time=linspace(t0,ts(end),10+2*K);
X=[X P(a1,time(K:end-K)-t0)];
Y=[Y P(a2,time(K:end-K)-t0)];
Ax=[Ax ddP(a1,time(K:end-K)-t0)];
Ay=[Ay ddP(a2,time(K:end-K)-t0)];
end
Max_a=max(am)
Mean_a=mean(am)
std_a=std(am)

%%
% close all
figure
scale=min(sqrt((X(2:end)-X(1:end-1)).^2+((Y(2:end)-Y(1:end-1)).^2)));
ColorFieldV(X,Y,Ax,Ay,scale,30)

figure
COL=jet(length(X)+1);
A=sqrt(Ax.^2+Ay.^2);

coli=floor((A-min(A))./(max(A)-min(A))*length(X))+1
COL(coli,:);

hold on
for J=1:length(X)
plot3(X(J),Y(J),A(J),'.','color',COL(coli(J),:),'markersize',15)
end


function plotmyfigures(ts,xs,vx,ax,ys,vy,ay,P,dP,ddP,a1,a2,t0,i)
figure
plot(ts,xs,'markersize',10)
hold on
plot(ts,P(a1,ts-t0),'linewidth',2)
title(['Difference to Linear Fit Over Time'])
xlabel('t [s]')
ylabel('\Deltax [m]')
set(gca,'fontsize',30)

figure
plot(ts,xs-P(a1,ts-t0),'markersize',10)
hold on
plot(ts,0,'linewidth',2)
title(['Difference to Linear Fit Over Time'])
xlabel('t [s]')
ylabel('\Delta\Deltax [m]')
set(gca,'fontsize',30)



figure
plot(ts,vx,'.','markersize',10)
hold on
plot(ts,dP(a1,ts-t0),'linewidth',2)
title(['Velocity Over Time'])
xlabel('t [s]')
ylabel('v [m]')
set(gca,'fontsize',30)

figure
plot(ts,ax,'.','markersize',10)
hold on
plot(ts,ddP(a1,ts-t0),'linewidth',2)
title('Acceleration Over Time')
xlabel('t [s]')
ylabel('a [m/s^2]')
set(gca,'fontsize',30)

figure
plot(ts,ys,'.')
hold on
plot(ts,P(a2,ts-t0))
title(['y (',num2str(i),')'])

figure
plot(ts,vy,'.')
hold on
plot(ts,dP(a2,ts-t0))
title(['vy (',num2str(i),')'])
set(gca,'fontsize',30)

figure
plot(ts,ay,'.')
hold on
plot(ts,ddP(a2,ts-t0))
title(['ay (',num2str(i),')'])
set(gca,'fontsize',30)
end
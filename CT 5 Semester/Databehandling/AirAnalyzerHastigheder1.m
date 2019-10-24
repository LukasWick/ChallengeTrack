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

% A=load('driftmåling1.txt');
A=load('data2_drag.txt');

t=A(:,1)/8;
x=A(:,2);
y=A(:,3);
% 
% A=load('driftmåling2.txt');
% t=[t;A(:,1)/8];
% x=[x;A(:,2)];
% y=[y;A(:,3)];

plot3(x,y,1:length(x))
figure 


% plot(x(1e4:end),'.')
% hold on
% plot(y(1e4:end),'.')


plot(x,'.')
hold on
plot(y,'.')


%%
% close all;
% i1=[3 820 1605 2395 3108 3871 4552 4782 5698 6205 7153 8801 9542 10865 11864 12076 12982 13187];
% i2=[210 957 1744 2539 3257 4005 4650 4923 5817 6324 7277 8897 9637 10972 11954 12199 13070 13320];
i1=[3 820 1605 2395 3108 3871 4552 4782 5698 6205 7153 8801 9542 10865 11864 12076 13187];
i2=[210 957 1744 2539 3257 4005 4650 4923 5817 6324 7277 8897 9637 10972 11954 12199 13320];

i1 = [4553,4653]
i2 = [4648,4779] 

% % i1 = [6037,6093]
% % i2 = [6087,6201] 
% % 
% % i1 = [6093,6209]
% % i2 = [6201,6335] 


% i1 = [10765,10865]
% i2 = [10855,10965]

% 
i1 = [7715,8045]
i2 = [8035,8380]
% 

% i1 = [6214,6349]
% i2 = [6329,6492]

% i1 = [12770,12920]
% i2 = [12910,13090]

am=[];
X=[];
Y=[];
Ax=[];
Ay=[];
for i=1:2
ts=t(i1(i):i2(i));

xs=x(i1(i):i2(i));
ys=y(i1(i):i2(i));
plot(xs,'.')
hold on
plot(ys,'.')

vx=gradient(xs,ts);
vy=gradient(ys,ts);
ax=gradient(vx,ts);
ay=gradient(vy,ts);

t0=ts(1);
a1(i,:)=nlinfit(ts,xs,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
a2(i,:)=nlinfit(ts,ys,@(a,t)P(a,t-t0),zeros(1,n_poly+1));
end
am=max(am)

t1=t(i1(1):i2(1))-t(i1(1));
for i = (i1(1):i2(1))-i1(1)+1
   t2(i) = fzero(@(t2) P(a1(1,:),t1(i))-P(a1(2,:),t2),t1(max(1,length(t1)-i)));     
end
t2 = t2';
Alpha = -(ddP(a1(2,:),t2)-ddP(a1(1,:),t1))./(dP(a1(2,:),t2)-dP(a1(1,:),t1));
figure
title('Alpha')
hold on
plot(Alpha)
figure
hold on
%%
plot(P(a1(2,:),t2),P(a2(2,:),t2))
plot(P(a1(1,:),t1),P(a2(1,:),t1))
plot(x(i1(1):i2(1)),y(i1(1):i2(1)),'.')
plot(x(i1(2):i2(2)),y(i1(2):i2(2)),'.')

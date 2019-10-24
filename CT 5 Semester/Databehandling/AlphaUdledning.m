clc; clear all; close all;
ns = [6,7,8,9,10,11]
alpha = []
for iindex = 1:6
n_poly = ns(iindex)
% n_poly=10;

Rfun=@(a,x)a(1);
dRfun=@(a,x) 0;
ddRfun=@(a,x) 0;
for i=1:n_poly
Rfun=@(a,x) Rfun(a,x)+a(i+1)*x.^i;
dRfun=@(a,x) dRfun(a,x)+i*a(i+1)*x.^(i-1);
end
for i=2:n_poly
ddRfun=@(a,x) ddRfun(a,x)+i*(i-1)*a(i+1)*x.^(i-2);
end



% A=load('driftmåling1.txt');
A=load('drag.txt');

t=A(:,1)/8;
x=A(:,2);
y=A(:,3);
dt = 1/240;
t = dt.*(1:length(t));
t = t';


% plot3(x,y,1:length(x))

% Define bounderys 
i1 = [1,63];
i2 = [61,139];

t1=t(i1(1):i2(1))-t(i1(1));
x1=x(i1(1):i2(1));
y1=y(i1(1):i2(1));

betaX1=nlinfit(t1,x1,@(a,t)Rfun(a,t1),ones(1,n_poly+1));
betaY1=nlinfit(t1,y1,@(a,t)Rfun(a,t1),ones(1,n_poly+1));

t2=t(i1(2):i2(2))-t(i1(2));
x2=x(i1(2):i2(2));
y2=y(i1(2):i2(2));

betaX2=nlinfit(t2,x2,@(a,t)Rfun(a,t2),ones(1,n_poly+1));
betaY2=nlinfit(t2,y2,@(a,t)Rfun(a,t2),ones(1,n_poly+1));


% for i = (i1(1):i2(1))-i1(1)+1
%    t2Adapted(i) = fzero(@(t2) Rfun(betaX1,t1(i))-Rfun(betaX2,t2),t1(max(1,length(t1)-i)));     
% end

timelengt = 1000;
ts1 = linspace(t1(1),t1(end),1000);
for i = 1:timelengt
   t2Adapted(i) = fzero(@(t2) Rfun(betaX1,ts1(i))-Rfun(betaX2,t2),ts1(max(1,length(ts1)-i)));     
end

Alpha = -(ddRfun(betaX2,t2Adapted)-ddRfun(betaX1,ts1))./(dRfun(betaX2,t2Adapted)-dRfun(betaX1,ts1));

%% x,y plot
% figure
% hold on 
% title(['Path on table '])
% xlabel('x [m]')
% ylabel('y [m]')
% set(gca,'fontsize',15)
% plot(x1,y1,'.')
% plot(x2,y2,'.')
% plot(Rfun(betaX1,ts1),Rfun(betaY1,ts1))
% % plot(Rfun(betaX2,ts2),Rfun(betaY2,ts2))
% plot(Rfun(betaX2,t2Adapted),Rfun(betaY2,t2Adapted))
% legend('Path 1','Path 2','Fit 1','Fit 2')

%% t,x plot
figure
hold on 
title(['Positions as function of time'])
xlabel('t [s]')
ylabel('x [m]')
set(gca,'fontsize',15)
plot(t1,x1,'.')
plot(t2,x2,'.')
plot(ts1,Rfun(betaX1,ts1))
plot(t2Adapted,Rfun(betaX2,t2Adapted))
legend('Path 1','Path 2','Fit 1','Fit 2')


%% x,v plot
% figure
% hold on 
% title(['Speed as function of position'])
% xlabel('x [m]')
% ylabel('v [m/s]')
% set(gca,'fontsize',15)
% plot(Rfun(betaX1,ts1),dRfun(betaX1,ts1))
% plot(Rfun(betaX2,t2Adapted),dRfun(betaX2,t2Adapted))
% legend('Fit 1','Fit 2')

%% x,a plot
figure
hold on 
title(['Acceleration as function of position'])
xlabel('x [m]')
ylabel('a [m/s^2]')
set(gca,'fontsize',15)
plot(Rfun(betaX1,ts1),ddRfun(betaX1,ts1))
plot(Rfun(betaX2,t2Adapted),ddRfun(betaX2,t2Adapted))
legend('Fit 1','Fit 2')

%% x,alpha plot
figure
hold on 
title(['Alpha as function of position'])
xlabel('x [m]')
ylabel('\alpha [1/s]')
set(gca,'fontsize',15)
% plot(Rfun(betaX1,ts1),Alpha)
xs= Rfun(betaX1,ts1);
plot(xs((timelengt/50):(timelengt-timelengt/50)),Alpha((timelengt/50):(timelengt-timelengt/50)))

%% alpha histogram
figure
hold on 
title(['Histogram of Alpha'])
xlabel('\alpha [1/s]')
set(gca,'fontsize',15)
histogram(Alpha((timelengt/50):(timelengt-timelengt/50)))
% %%
% 
% figure
% T = t(50)
% tList = ones([1,n_poly+1])
% for i = 2:n_poly+1
%     tList = [tList((1:i-1)),tList((i:end)).*T];
% end
% plot(betaX1.*tList)
alpha(iindex,:)=Alpha;
end
figure
hold on 
title(['Alpha as function of position'])
xlabel('x [m]')
ylabel('\alpha [1/s]')
set(gca,'fontsize',15)
% plot(Rfun(betaX1,ts1),Alpha)
plot(alpha')

%% alpha histogram
figure
hold on 
title(['Histogram of Alpha'])
xlabel('\alpha [1/s]')
set(gca,'fontsize',15)
histogram(alpha)

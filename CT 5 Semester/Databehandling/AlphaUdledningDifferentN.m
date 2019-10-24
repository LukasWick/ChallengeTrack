clc; clear all; close all;

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

dt = 1/240
t1= (t(i1(1):i2(1))-t(i1(1)));
t1 = dt.*(1:length(t1));
t1 = t1';
x1=x(i1(1):i2(1));
y1=y(i1(1):i2(1));

t2=t(i1(2):i2(2))-t(i1(2));
t2 = (dt.*(1:length(t2)));
t2 = t2';
x2=x(i1(2):i2(2));
y2=y(i1(2):i2(2));

ns = [6,7,8,9,10,11]
alpha = []
for iindex = 1:6
n_poly = ns(iindex)


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

linFun = @(a,t) a(1)+a(2).*t

linBetaX1 = nlinfit(t1,x1,linFun,[1,1])
linBetaX2 = nlinfit(t2,x2,linFun,[1,1])

deltaX1 = x1 - linFun(linBetaX1,t1)
deltaX2 = x2 - linFun(linBetaX2,t2)

betaX1=nlinfit(t1,deltaX1,@(a,t)Rfun(a,t),ones(1,n_poly+1));
betaX2=nlinfit(t2,deltaX2,@(a,t)Rfun(a,t),ones(1,n_poly+1));

x1Fit =@(t) linFun(linBetaX1,t)+Rfun(betaX1,t)
x2Fit =@(t) linFun(linBetaX2,t)+Rfun(betaX2,t)
v1Fit =@(t) linBetaX1(2)+dRfun(betaX1,t)
v2Fit =@(t) linBetaX2(2)+dRfun(betaX2,t)

timelengt = 1000;
ts1 = linspace(t1(1),t1(end),1000);
for i = 1:timelengt
   ts2(i) = fzero(@(t2) x1Fit(ts1(i))-x2Fit(t2),ts1(max(1,length(ts1)-i)));     
end

Alpha = -(ddRfun(betaX2,ts2)-ddRfun(betaX1,ts1))./(v1Fit(ts1)-v2Fit(ts2));

% x,y plot
figure
hold on 
title(['Path on table '])
xlabel('x [m]')
ylabel('y [m]')
set(gca,'fontsize',15)
plot(x1,y1,'.')
plot(x2,y2,'.')
plot(x1Fit(t1),y1)
plot(x2Fit(t2),y2)

% plot(Rfun(betaX1,ts1),Rfun(betaY1,ts1))
% % plot(Rfun(betaX2,ts2),Rfun(betaY2,ts2))
% plot(Rfun(betaX2,t2Adapted),Rfun(betaY2,t2Adapted))
legend('Path 1','Path 2','Fit 1','Fit 2')

%% t,deltax plot
figure
hold on 
title(['Residuel as function of time'])
xlabel('t [s]')
ylabel('\Deltax [m]')
set(gca,'fontsize',15)
plot(t1,x1-linFun(linBetaX1,t1),'.')
plot(t2,x2-linFun(linBetaX2,t2),'.')
plot(ts1,Rfun(betaX1,ts1))
plot(ts2,Rfun(betaX2,ts2))
legend('Path 1','Path 2','Fit 1','Fit 2')


%% t,x plot
figure
hold on 
title(['Positions as function of time'])
xlabel('t [s]')
ylabel('x [m]')
set(gca,'fontsize',15)
plot(t1,x1,'.')
plot(t2,x2,'.')
plot(ts1,x1Fit(ts1))
plot(ts2,x2Fit(ts2))
legend('Path 1','Path 2','Fit 1','Fit 2')


%% x,v plot
figure
hold on 
title(['Speed as function of position'])
xlabel('x [m]')
ylabel('v [m/s]')
set(gca,'fontsize',15)
plot(x1Fit(ts1),v1Fit(ts1))
plot(x2Fit(ts2),v2Fit(ts2))
legend('Fit 1','Fit 2')

%% x,a plot
figure
hold on 
title(['Acceleration as function of position'])
xlabel('x [m]')
ylabel('a [m/s^2]')
set(gca,'fontsize',15)
plot(x1Fit(ts1),ddRfun(betaX1,ts1))
plot(x2Fit(ts2),ddRfun(betaX2,ts2))
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
plot(x1Fit(ts1),Alpha)

%% alpha histogram
figure
hold on 
title(['Histogram of Alpha'])
xlabel('\alpha [1/s]')
set(gca,'fontsize',15)
histogram(Alpha,100)
[N,E] = histcounts(Alpha,100)
edges = E(1:(end-1))+diff(E)/2
edges = [edges,edges+edges(end)-edges(1)]
N = [N,zeros(size(N))]

edges = edges(50:(end-50))
N = N(50:(end-50))
% plot(edges,N,'.')
gauss = @(beta,x) beta(1).*exp(-((x-beta(2))./beta(3)).^2./2)
beta = nlinfit(edges,N,gauss,[max(N),0,2])
plot(edges,gauss(beta,edges))
legend('Alpha his',['Gausfit \alpha = ',num2str(beta(2))])

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
histogram(alpha,100)
[N,E] = histcounts(alpha,100)
edges = E(1:(end-1))+diff(E)/2
edges = [edges,edges+edges(end)-edges(1)]
N = [N,zeros(size(N))]

edges = edges(50:(end-50))
N = N(50:(end-50))
% plot(edges,N,'.')
gauss = @(beta,x) beta(1).*exp(-((x-beta(2))./beta(3)).^2./2)
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(edges,N,gauss,[max(N),0,2])
plot(edges,gauss(beta,edges))
legend('Alpha his',['Gausfit \alpha = ',num2str(beta(2))])

ci = nlparci(beta,R,'jacobian',J)

alp = beta(2)
alpUs = norm((ci(2,1)-ci(2,2)))/2

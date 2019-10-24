clc; clear all; close all;
%%
alphaSet = 0.1
stoj = 1e-3
T = 0.2500;
n = 60;
v0 = 2.3
accWidth = 3
acc_n = 5
ns = 6:11;
number_ofmeasurements = 1


alpha = []
for LongIndex = 1:number_ofmeasurements
times = linspace(0,T,n);
dat0 = [0;v0];

acc_points = normrnd(0,accWidth,[acc_n,1]);
acc_points_x =linspace(0,T*dat0(2),acc_n);%[-0.2;-0.1;rand([3,1])*T*dat0(2)/2+T*dat0(2)/2,;T*dat0(2)+[0.1;0.2]];
accFun = @(x) spline(acc_points_x,acc_points,x);


[t,dat] = ode45(@(t,dat) diffFun(t,dat,accFun,alphaSet),times,dat0);
xs = dat(:,1);
vs = dat(:,2);

%%
% figure
% plot(xs,accFun(xs))
% hold on 
% plot(acc_points_x,acc_points,'.')


dat0 = [xs(end);-vs(end)];

[temp,dat] = ode45(@(t,dat) diffFun(t,dat,accFun,alphaSet),times,dat0);
t = [t;t(end)+temp];
xs = [xs;dat(:,1)];
vs = [vs;dat(:,2)];


%%
t=t;
x=xs+normrnd(0,stoj,size(t));
y=zeros(size(t))+normrnd(0,stoj,size(t));

% plot3(x,y,1:length(x))
% Define bounderys 
i1 = [1,61];
i2 = [59,120];

dt = 1/240;
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

for iindex = 1:length(ns)
n_poly = ns(iindex);


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

linFun = @(a,t) a(1)+a(2).*t;

linBetaX1 = nlinfit(t1,x1,linFun,[1,1]);
linBetaX2 = nlinfit(t2,x2,linFun,[1,1]);

deltaX1 = x1 - linFun(linBetaX1,t1);
deltaX2 = x2 - linFun(linBetaX2,t2);

betaX1=nlinfit(t1,deltaX1,@(a,t)Rfun(a,t),ones(1,n_poly+1));
betaX2=nlinfit(t2,deltaX2,@(a,t)Rfun(a,t),ones(1,n_poly+1));

x1Fit =@(t) linFun(linBetaX1,t)+Rfun(betaX1,t);
x2Fit =@(t) linFun(linBetaX2,t)+Rfun(betaX2,t);
v1Fit =@(t) linBetaX1(2)+dRfun(betaX1,t);
v2Fit =@(t) linBetaX2(2)+dRfun(betaX2,t);

timelengt = 1000;
ts1 = linspace(t1(1),t1(end),1000);
for i = 1:timelengt
   ts2(i) = fzero(@(t2) x1Fit(ts1(i))-x2Fit(t2),ts1(max(1,length(ts1)-i)));     
end

Alpha = -(ddRfun(betaX2,ts2)-ddRfun(betaX1,ts1))./(v1Fit(ts1)-v2Fit(ts2));

%% x,y plot
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
plot(x1Fit(ts1),accFun(x1Fit(ts1)))

legend('Fit 1','Fit 2','Model')

%% x,alpha plot
lenA = length(Alpha);
ALPHA = Alpha(lenA/10:(end-lenA/10));
xs = x1Fit(ts1);
XS = xs(lenA/10:(end-lenA/10));
figure
hold on 
title(['Alpha as function of position'])
xlabel('x [m]')
ylabel('\alpha [1/s]')
set(gca,'fontsize',15)
% plot(Rfun(betaX1,ts1),Alpha)
plot(XS,ALPHA)

%% alpha histogram
figure
hold on 
title(['Histogram of Alpha'])
xlabel('\alpha [1/s]')
set(gca,'fontsize',15)
histogram(ALPHA,100)
[N,E] = histcounts(ALPHA,100)
edges = E(1:(end-1))+diff(E)/2
plot(edges,N,'.')
gauss = @(beta,x) beta(1).*exp(-((x-beta(2))./beta(3)).^2./2)
beta = nlinfit(edges,N,gauss,[max(N),0,2])
plot(edges,gauss(beta,edges))
legend('Alpha his',['Gausfit \alpha = ',num2str(beta(2))])

 alpha(end+1,:)=ALPHA;
end
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
h=histogram(alpha,1000,'Normalization','probability');
h.EdgeColor = 'none';

[N,E] = histcounts(alpha,1000,'Normalization','probability');
edges = E(1:(end-1))+diff(E)/2;
dE = diff(E);

alphaFound = mean(alpha(:))
alphaFoundStd = std(alpha(:))


egess = linspace(edges(1),edges(end),1000);
% gauss = dE(1)*1/sqrt(2*pi*alphaFoundStd^2).*exp(-((egess-alphaFound)./alphaFoundStd).^2./2);

% plot(egess,gauss)
gauss = @(beta,x) beta(1).*exp(-((x-beta(2))./beta(3)).^2./2);
plot(edges,gauss([dE(1)*1/sqrt(2*pi*alphaFoundStd^2),alphaFound,alphaFoundStd],edges))

[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(edges,N,gauss,[dE(1)*1/sqrt(2*pi*alphaFoundStd^2),alphaFound,alphaFoundStd])
plot(edges,gauss(beta,edges))
% plot(edges,N)
legend('Alpha his',['mean \alpha = ',num2str(alphaFound)],['gausfit \alpha = ',num2str(beta(2))])

ci = nlparci(beta,R,'jacobian',J);
% 
alp = beta(2)
alpUs = norm((ci(2,1)-ci(2,2)))/2
% 

%%

function dDat = diffFun(t,dat,accFun,alpha)
x = dat(1);
v = dat(2);
dxdt = v;
dvdt = accFun(x)-alpha*v;
dDat = [dxdt;dvdt];
end


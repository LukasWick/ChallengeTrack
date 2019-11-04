clc; clear all; close all;
%%
alphaSet = 0.1;
stoj = 2e-4;

L = 0.6;
fps = 240;
accWidth = 0.5
changeInAcc = 0.1
acc_n = 4
ns = 5:9;
number_ofmeasurements = 20;
v0_const = 1
% T = 0.2500;


alpha = [];
%Genner punkter accelerationskurven følger
acc_points_const = normrnd(0,accWidth,[acc_n,1]);
acc_points_x = linspace(0,L,acc_n);

vxs=[];
vys=[];
axs=[];
ays=[];


for i = 1:number_ofmeasurements

v0 = normrnd(0,v0_const);
v0 = v0 + sign(v0)*(abs(normrnd(0,v0_const/10))+0.5);
if (v0>0)
    dat0 = [0;v0]; %x0;v0
else
    dat0 = [L;v0]; %x0;v0
end
T = L/abs(v0);
n = round(T*fps);
times = linspace(0,T,n);
if (v0>0)
    dat0 = [0;v0]; %x0;v0
else
    dat0 = [L;v0]; %x0;v0
end

% Ændre en lille smule på accelerationen (mere jo længere fra stød kanten)
acc_points =acc_points_const+ normrnd(0,accWidth*changeInAcc,[acc_n,1]);
accFun = @(x) spline(acc_points_x,acc_points,x);

%Løs diff ud
    [t,dat] = ode45(@(t,dat) diffFun(t,dat,accFun,alphaSet),times,dat0);
    Xs = dat(:,1);
    vs = dat(:,2);
    t=t;
    x=Xs+normrnd(0,stoj,size(t));
    y=[t*v0/10]+normrnd(0,stoj,size(t));
    


    vx=gradient(x,t);
    vy=gradient(y,t);
    ax=gradient(vx,t);
    ay=gradient(vy,t);
    
    vxs=[vxs; vx];
    vys=[vys; vy];
    axs=[axs; ax];
    ays=[ays; ay];
end
plot(vxs,axs,'.')
xlabel('vx [m/s]')
ylabel('ax [m/s^2]')

% vs=sqrt(vxs.^2+vys.^2);
% as=sqrt(axs.^2+ays.^2);
% plot(vs,as,'.')

[drags,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(vxs,axs,@(drags,v) -drags(1)*v+drags(2),[0.1 1]);

alphax=drags(1)
ci = nlparci(drags,R,'jacobian',J);
usx = (ci(1,2)-ci(1,1))/4


figure
plot(vys,ays,'.')
xlabel('vy [m/s]')
ylabel('ay [m/s^2]')


[drags,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(vys,ays,@(drags,v) -drags(1)*v+drags(2),[0.1 1]);

alphay=drags(1)
ci = nlparci(drags,R,'jacobian',J);
usy = (ci(1,2)-ci(1,1))/4


alpha=(alphax/usx^2+alphay/usy^2)/(1/usx^2+1/usy^2)
us=1/sqrt(1/usx^2+1/usy^2)

mean_v=mean(sqrt(vxs.^2+vys).^2)
mean_drag=alpha*mean_v


function dDat = diffFun(t,dat,accFun,alpha)
x = dat(1);
v = dat(2);
dxdt = v;
dvdt = accFun(x)-alpha*v;
dDat = [dxdt;dvdt];
end


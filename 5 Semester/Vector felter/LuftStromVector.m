clear all; close all;clc;
addpath('robustDifferentiators')
% addpath(genpath('robustDifferentiators'))
delimiter = ';';
startRow = 3;
formatSpec = '%f%f%f%[^\n\r]';
filenames = {'data.txt','data2_drag.txt'};
t = [];
x = [];
y = [];
normaliseringsfactor =[1,8]
figure
for i = 1:2
    fileID = fopen(filenames{i},'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    T = dataArray{:, 1};
    X = dataArray{:, 2};
    Y = dataArray{:, 3};
    T = T(1000:(end-1000))*normaliseringsfactor(i);
    X = X(1000:(end-1000));
    Y = Y(1000:(end-1000));
    t = [t;T];
    x = [x;X];
    y = [y;Y];
    plot(T)
    hold on
end
% clearvars filename delimiter startRow formatSpec fileID dataArray ans X Y T;

figure


plot(x,y)
axis equal


xmax = 0.4
xmin = 0.05
ymax = 0.3
ymin = 0.06

% xmax = max(x)
% xmin = min(x)
% ymax = max(y)
% ymin = min(y)


dxdt = gradient(x,t);
d2xdt2 = gradient(dxdt,t);
dydt = gradient(y,t);
d2ydt2 = gradient(dydt,t);

v = sqrt(dxdt.^2+dydt.^2)
a = sqrt(d2xdt2.^2+d2ydt2.^2)

acutoff = 0.6
vcutoff = 0.5
t = t((v<vcutoff)&(a<acutoff))
x = x((v<vcutoff)&(a<acutoff))
y = y((v<vcutoff)&(a<acutoff))
dxdt = dxdt((v<vcutoff)&(a<acutoff))
d2xdt2 = d2xdt2((v<vcutoff)&(a<acutoff))
dydt = dydt((vcutoff)&(a<acutoff))
d2ydt2 = d2ydt2((v<vcutoff)&(a<acutoff))

figure 
plot(x,y)
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])


N=6
dx = (xmax-xmin)/N
dy = (ymax-ymin)/N
luftmodstand =0
figure
for i=1:N
    x0 = xmin+(i-1)*dx;
    x1 = xmin+(i)*dx;
    columIndex = ((x0<x)&(x<x1));
    X_colum = x(columIndex);
    Y_colum = y(columIndex);
    d2xdt2_colum = d2xdt2(columIndex);
    d2ydt2_colum = d2ydt2(columIndex);
    dxdt_colum = dxdt(columIndex);
    dydt_colum = dydt(columIndex);
    t_colum = t(columIndex);
    for j =1:N
        y0 = ymin+(j-1)*dy;
        y1 = ymin+(j)*dy;
        rowIndex = (y0<Y_colum)&(Y_colum<y1);
        X_row = X_colum(rowIndex);
        Y_row = Y_colum(rowIndex);
        d2xdt2_row = d2xdt2_colum(rowIndex);
        d2ydt2_row = d2ydt2_colum(rowIndex);
        dxdt_row = dxdt_colum(rowIndex);
        dydt_row = dydt_colum(rowIndex);
        t_row = t_colum(rowIndex);
        hold on
        plot(t_row,d2xdt2_row,'.')
        fieldX(j,i)=mean(d2xdt2_row-dxdt_row*luftmodstand);
        fieldY(j,i)=mean(d2ydt2_row-dydt_row*luftmodstand);
    end
end
 Xs = linspace(xmin+dx,xmax-dx,N);
 Ys = linspace(ymin+dy,ymax-dy,N);
% 
figure
[xs,ys] = meshgrid(Xs,Ys);
quiver(xs,ys,fieldX,fieldY)
%%

% alphaX = fzero(@(alpha)mean((d2xdt2-alpha.*dxdt).*dxdt),0.1)
% alphaY = fzero(@(alpha)mean((d2ydt2-alpha.*dydt).*dydt),0.1)
% 
% alpha = linspace(-0.1,0.1,1000)
% figure
% plot(alpha,mean((d2xdt2-alpha.*dxdt).*dxdt))
% hold on
% plot(alpha,mean((d2ydt2-alpha.*dydt).*dydt))
% v = sqrt(dxdt.^2+dydt.^2)

% alphaV = fminsearch(@(alpha) mean(((d2xdt2-alpha.*dxdt).*dxdt).^2+((d2ydt2-alpha.*dydt).*dydt).^2),0.1)

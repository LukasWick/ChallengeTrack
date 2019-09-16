clear all; close all;clc;
addpath('robustDifferentiators')
% addpath(genpath('robustDifferentiators'))
delimiter = ';';
startRow = 3;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen('data2_drag.txt','r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
t = dataArray{:, 1};
x = dataArray{:, 2};
y = dataArray{:, 3};
t = t(1:(end-0)).*8
x = x(1:(end-0));
y = y(1:(end-0));
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

plot(x,y)
axis equal

xmax = 0.4
xmin = -0.1
ymax = 0.4
ymin = 0.2

figure 
plot(x,y)
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])

dxdt = gradient(x,t);
d2xdt2 = gradient(dxdt,t);
dydt = gradient(y,t);
d2ydt2 = gradient(dydt,t);

figure
v = sqrt(dxdt.^2+dydt.^2)
a = sqrt(d2xdt2.^2+d2ydt2.^2)
plot(v)
figure
v =  [v(2750:5734);v(6093:8674)]
t =  [t(2750:5734);t(6093:8674)]
a =  [a(2750:5734);a(6093:8674)]
plot(v)
figure
plot(a)
title('a')
%%

prevDelim = 1
veloseties ={}
times = {}
for i = 2:length(v)
   if(a(i)>0.35|(t(i)-t(i-1))>1)
       if prevDelim ~=(i-1)
            veloseties = cat(1,veloseties,v(prevDelim:i-1))
            times = cat(1,times,t(prevDelim:i-1))
            prevDelim = i+1
       end
   end
end
veloseties = cat(1,veloseties,v(prevDelim:length(v)))
times = cat(1,times,t(prevDelim:length(v)))

figure
hold on

plot(t,v)
[n,m] =size(times)
for i = 1:n
    beta(:,i) = nlinfit(times{i},veloseties{i},@(beta,t) beta(1).*t+beta(2),[1,0]);
    alpha(i) = beta(1,i)./mean(veloseties{i})
    plot(times{i},veloseties{i},'*')
end
figure
plot(beta(1,:))

mean(alpha)
std(alpha)/sqrt(n)


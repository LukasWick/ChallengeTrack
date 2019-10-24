% clear all; close all;clc;
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
t = t(1:(end-0));
x = x(1:(end-0));
y = y(1:(end-0));
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
figure
plot(t,x,'.')

figure
indexes = (t>175.5/8&t<185.4/8) 
x = x(indexes)
t = t(indexes)
y = y(indexes)
plot(t,x,'.')

figure
beta = nlinfit(t,x,@(beta,t) beta(1).*t+beta(2),[1,0])
plot(t,x-beta(1).*t-beta(2),'.')

figure
plot(x,y)
%%
figure
plot(x,y)
axis equal

xmax = 0.5
xmin = -0.2
ymax = 0.5
ymin = 0.3

xmax = 2
xmin = -2
ymax = 2
ymin = -2


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
v = sqrt(dxdt.^2+dydt.^2);
a = sqrt(d2xdt2.^2+d2ydt2.^2);
A = gradient(v,t);

plot(v)
figure
v =  [v(2750:5734);v(6093:8674)];
t =  [t(2750:5734);t(6093:8674)];
a =  [a(2750:5734);a(6093:8674)];
A =  [A(2750:5734);A(6093:8674)]
plot(v)
figure
plot(a)
title('a')
figure
plot(A)

title('A')
figure

%%

prevDelim = 1
veloseties ={}
times = {} 
As = {}
for i = 2:length(v)
   if(a(i)>0.35|(t(i)-t(i-1))>1|abs(A(i))>0.35| x(i)>xmax |x(i)<xmin|y(i)>ymax |y(i)<ymin)
       if prevDelim ~=(i-1)
            veloseties = cat(1,veloseties,v(prevDelim:i-1))
            As = cat(1,As,A(prevDelim:i-1))
            times = cat(1,times,t(prevDelim:i-1))
            prevDelim = i+1
       end
   end
end
veloseties = cat(1,veloseties,v(prevDelim:length(v)))
times = cat(1,times,t(prevDelim:length(v)))
As = cat(1,As,A(prevDelim:length(v)))

figure(9)
hold on

plot(t,v)
[n,m] =size(times)
Ases = []
    figure(9)
beta = []
for i = 1:n

    beta(:,i) = nlinfit(times{i},log(veloseties{i}),@(beta,t) beta(1).*t+beta(2),[1,0]);
    alpha(i) = beta(1,i)%./mean(veloseties{i});
    plot(times{i},log(veloseties{i}),'*');
end
figure(10) 

for i = 1:n
hold on
plot(times{i},As{i},'*')
Ases = [Ases;As{i}]
end
figure
plot(beta(1,:))

mean(alpha)
std(alpha)/sqrt(n)


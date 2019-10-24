function ColorFieldV(x,y,fx,fy,Scale)
a=length(x);
Scale;
d=sqrt(fx.^2+fy.^2);
L=max(d);
g=min(d);
e=L-g;
e=e+(e==0);

subplot(1,6,1)
hold on
N=500;
a1=jet(N);
for i=1:N
fill([-0.1 -0.1 0.1 0.1],[L/N*(i-1) L/N*i L/N*i L/N*(i-1)],a1(i,:),'linestyle','none')
end
axis([-0.1 0.1 0 L])
xlabel('')
ylabel('Length of vector')
title('Color Index')
set(gca,'xtick',[])

subplot(1,6,[2 3 4 5 6])
color=jet(a);
for j=1:a  
D=d(j)/Scale;
f=(D*Scale-g)/e;
f=round(f*a);
f=f+(f==0);
arrow(x(j),y(j),fx(j)./D,fy(j)./D,f,Scale,color)
end
axis equal
title('Vector Field')
xlabel('x')
ylabel('y')

end
function arrow=arrow(x,y,fx,fy,f,Scale,color)
h=0.2;
theta = 0.5;
R1 = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
R2 = [cos(theta)  sin(theta) ; -sin(theta)  cos(theta)]; 
xe=x+fx;
ye=y+fy;
v1 = h*([x y]-[xe ye]);
v2 = v1*R1;
v3 = v1*R2;
x1 = [xe ye] + v2;
x2 = [xe ye] + v3;
plot([xe x2(1)],[ye x2(2)],'LineWidth',Scale,'Color',color(f,:))
hold on;
plot([x xe x1(1)],[y ye x1(2)],'LineWidth',Scale,'Color',color(f,:));
end
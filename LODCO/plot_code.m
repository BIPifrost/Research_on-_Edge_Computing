x = [10,19];
plot(x);

t = 0:0.02:2*pi;
y1 = sin(t);
y2 = sin(2*t);
y3 = sin(0.5*t);
plot(t,y1,t,y2,t,y3);

x = [2021,2022,2023];
y = [10,20,30];
bar(x,y);
x = [2021,2022,2023];
y = [10,20;30,40;50,60];
bar(x,y);

x=randn(1000,1);
nbins = 25;
histogram(x,nbins);

x = [3,2,5,10];
pie(x);

x=[1 3 5 6  10 20]
y=[ 13 45 73 23 78 23]
scatter(x,y,'red','filled')

t = 0:0.02:2*pi;
y1 = sin(t);
y2 = sin(2*t);
y3 = sin(0.5*t);
plot(t,y1,t,y2,t,y3);
axis([0,6.5,-1.5,1.5])
title('Three Sine Functions','Fontsize',12);
xlabel('X');
ylabel('Y');
legend('sin(x)','sin(2x)','sin(0.5x)')

text(2.5, sin(2.5),'sin(x)');
text(2.5, sin(5),'sin(2x)');

figure('Name','Subplot Demo');
subplot(2,2,1);
x=randn(1000,1);
nbins = 25;
histogram(x,nbins);

subplot(2,2,2);
t = 0:0.02:2*pi;
y1 = sin(t);
y2 = sin(2*t);
y3 = sin(0.5*t);
plot(t,y1,t,y2,t,y3);

subplot(2,2,3);
x = [3,2,5,10];
pie(x);

subplot(2,2,4);
x = [10,19];
plot(x);

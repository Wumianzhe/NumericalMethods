dir = '~/Documents/Numerical/build/lab5/';
file = fopen([dir,'vals.bin'],'rb');
data = fread(file,[4,Inf],'double');
fclose(file);

a = pi / 4;
b = 2 * pi;
n = 10;

figure;
hold on;
grid on;
title('Графики функции и приближений');
axis([a-0.05,b+0.05,-8,2]);
plot(data(1,:),data(2,:),"b;func;");
plot(data(1,:),data(3,:),"k;y(H);");
plot(data(1,:),data(4,:),"r;y(H/2);");

figure;
hold on;
grid on;
axis([a-0.05,b+0.05]);
title('Ошибка приближений');
plot(data(1,:),data(3,:)-data(2,:),"k;\Deltay(H);");
plot(data(1,:),data(4,:)-data(2,:),"r;\Deltay(H/2);");

file = fopen([dir,'errs.bin'],'rb');
errdata = fread(file,[15,Inf],'double');
fclose(file);

figure;
hold on;
grid on;
axis([a-0.05, b+0.05])
title('Все ошибки');
for i = 1:15
  ax = linspace(a,b,n+1);
  semilogy(ax, errdata(i,:), "color", 1./15.*[i,0,i]);
endfor

data = importdata([dir, "errs.csv"],',')';
figure;
hold on;
grid on;
h = (b-a)./(n*2.^data(2,:));
title('Ошибка от заданной точности');
loglog(data(1,:),data(3,:),"k;err;");
loglog(data(1,:),data(1,:),"b;eps;");
loglog(data(1,:),errdata(:,2),"r;eps[1];");

figure;
hold on;
grid on;
title('Итерации от заданной точности');
semilogx(data(1,:),data(2,:),"k;iter;");

data = importdata([dir, "delta.csv"],',')';
figure;
title('Ошибка от возмущения входных данных');
loglog(data(1,:),data(2,:),"k;err;");

data = importdata([dir, "step.csv"],',')';
figure;
title('Ошибка от длины шага');
hold on;
grid on;
loglog(data(1,:), data(2,:), "b;err;");
loglog(data(1,:), data(1,:).^4, "k;h^4;");


dir = '~/Documents/Numerical/build/lab7/';
file = fopen([dir,'vals.bin'],'rb');
data = fread(file,[4,Inf],'double');
fclose(file);

a = 1;
b = 2;

figure;
hold on;
grid on;
title('Графики функции и приближений');
axis([a-0.05,b+0.05]);
plot(data(1,:),data(2,:),"b;func;");
plot(data(1,:),data(3,:),"r;y(H);");
plot(data(1,:),data(4,:),"k;y(H/2);");


figure;
hold on;
grid on;
title('Графики ошибок на заданном отрезке');
plot(data(1,:),data(3,:) - data(2,:),"r;\Delta(H);");
plot(data(1,:),data(4,:) - data(2,:),"k;\Delta(H/2);");

data = importdata([dir, "errs.csv"],',')';

figure;
hold on;
grid on;
title('График ошибки от величины шага')
loglog(data(1,:),data(2,:),"k;err;");
loglog(data(1,:),data(1,:).^2,"b;h^2;")
coefs = data(2,:) ./ (data(1,:).^2);

data = importdata([dir, "delta.csv"],',')';

figure;
hold on;
grid on;
title('График ошибки от возмущения входных данных')
loglog(data(1,:),data(2,:),"r;+\Delta;");
loglog(data(1,:),data(3,:),"b;-\Delta;");

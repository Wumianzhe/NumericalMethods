file = fopen('~/Documents/Numerical/build/lab2/smooth.bin','rb');
data = fread(file,[5,Inf],'double');
fclose(file);

a = 1;
b = 2.5;
n_low = 7;
n_high = 1000;

figure;
hold on;
grid on;
title('Графики функции и приближений');
axis([a-0.05,b+0.05 1.5 5]);
plot(data(1,:),data(5,:),"b;func;",'LineWidth',1.25);
plot(data(1,:),data(2,:),"k;S_3;",'LineWidth',1.25);
plot(data(1,:),data(3,:),"g;S_7;",'LineWidth',1.25);
plot(data(1,:),data(4,:),"r;S_{10};",'LineWidth',1.25);

figure;
hold on;
grid on;
title('Ошибка для приближения функции');
axis([a,b]);
plot(data(1,:),data(2,:)-data(5,:),"k;\\DeltaS_3;",'LineWidth',1.25);
plot(data(1,:),data(3,:)-data(5,:),"g;\\DeltaS_7;",'LineWidth',1.25);
plot(data(1,:),data(4,:)-data(5,:),"b;\\DeltaS_{10};",'LineWidth',1.25);


file = fopen('~/Documents/Numerical/build/lab2/nsmooth.bin','rb');
data = fread(file,Inf,'double');
fclose(file);

figure;
hold on;
grid on;
title('Ошибка в зависимости от числа точек');
semilogy(n_low:n_high,data,'LineWidth',1.25);

file = fopen('~/Documents/Numerical/build/lab2/dsmooth.bin','rb');
data = fread(file, [2,Inf], 'double');
fclose(file);

figure;
hold on;
grid on;
axis([-1.05,2.05]);
title('Ошибка в зависимости от граничных условий');
semilogy(data(1,:),data(2,:),'LineWidth',1.25);




file = fopen('~/Documents/Numerical/build/lab2/corner.bin','rb');
data = fread(file,[5,Inf],'double');
fclose(file);

figure;
hold on;
grid on;
title('Графики негладкой функции и приближений');
axis([a-0.05,b+0.05 1.5 5]);
plot(data(1,:),data(5,:),"b;func;",'LineWidth',1.25);
plot(data(1,:),data(2,:),"k;S_3;",'LineWidth',1.25);
plot(data(1,:),data(3,:),"g;S_7;",'LineWidth',1.25);
plot(data(1,:),data(4,:),"r;S_{10};",'LineWidth',1.25);

figure;
hold on;
grid on;
title('Ошибка для приближения негладкой функции');
axis([a,b]);
plot(data(1,:),data(2,:)-data(5,:),"k;\\DeltaS_3;",'LineWidth',1.25);
plot(data(1,:),data(3,:)-data(5,:),"g;\\DeltaS_7;",'LineWidth',1.25);
plot(data(1,:),data(4,:)-data(5,:),"b;\\DeltaS_{10};",'LineWidth',1.25);


file = fopen('~/Documents/Numerical/build/lab2/ncorner.bin','rb');
data = fread(file,Inf,'double');
fclose(file);

figure;
hold on;
grid on;
title('Ошибка в зависимости от числа точек для негладкой ф-и');
semilogy(n_low:n_high,data,'LineWidth',1.25);

file = fopen('~/Documents/Numerical/build/lab2/dcorner.bin','rb');
data = fread(file, [2,Inf], 'double');
fclose(file);

figure;
hold on;
grid on;
axis([-1.05,2.05]);
title('Ошибка в зависимости от граничных условий для негладкой ф-и');
semilogy(data(1,:),data(2,:),'LineWidth',1.25);

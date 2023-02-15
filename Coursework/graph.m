dir = '~/Documents/Numerical/build/course/';
file = fopen([dir,'fpoly.bin'],'rb');
data_p = fread(file,[5,Inf],'double');
fclose(file);
file = fopen([dir,'fspline.bin'],'rb');
data_s = fread(file,[5,Inf],'double');
fclose(file);

a = 1;
b = 2.5;

figure;
subplot(1,2,1);
hold on;
grid on;
title('Графики функции и приближений полиномами');
axis([a-0.05,b+0.05 1.5 5]);
plot(data_p(1, :),data_p(5,:),"b;func;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(2,:),"k;P_3;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(3,:),"g;P_5;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(4,:),"r;P_7;", 'LineWidth',1.25);
legend(legend("show"), "location", "northwest");

subplot(1,2,2);
hold on;
grid on;
title('Графики функции и приближений сплайнами');
axis([a-0.05,b+0.05 1.5 5]);
plot(data_s(1, :),data_s(5,:),"b;func;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(2,:),"k;S_3;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(3,:),"g;S_5;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(4,:),"r;S_7;", 'LineWidth',1.25);
legend(legend("show"), "location", "northwest");

figure;
subplot(1,2,1);
hold on;
grid on;
title('Графики ошибки приближений полиномами');
axis([a-0.05,b+0.05]);
plot(data_p(1, :),data_p(2,:)-data_p(5,:),"k;P_3;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(3,:)-data_p(5,:),"g;P_5;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(4,:)-data_p(5,:),"r;P_7;", 'LineWidth',1.25);

subplot(1,2,2);
hold on;
grid on;
title('Графики ошибки приближений сплайнами');
axis([a-0.05,b+0.05]);
plot(data_s(1, :),data_s(2,:)-data_p(5,:),"k;S_3;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(3,:)-data_p(5,:),"g;S_5;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(4,:)-data_p(5,:),"r;S_7;", 'LineWidth',1.25);

file = fopen([dir,'nfpoly.bin'],'rb');
data_p = fread(file,Inf,'double');
fclose(file);
file = fopen([dir,'nfspline.bin'],'rb');
data_s = fread(file,Inf,'double');
fclose(file);

figure;
hold on;
grid on;
title('Зависимость ошибки от числа точек для гладкой функции');
semilogy(3:size(data_p)(1)+2,data_p,'k;poly;');
semilogy(3:size(data_s)(1)+2,data_s,'b;spline;');

file = fopen([dir,'gpoly.bin'],'rb');
data_p = fread(file,[5,Inf],'double');
fclose(file);
file = fopen([dir,'gspline.bin'],'rb');
data_s = fread(file,[5,Inf],'double');
fclose(file);

a = 1;
b = 2.5;

figure;
subplot(1,2,1);
hold on;
grid on;
title('Графики функции и приближений полиномами');
axis([a-0.05,b+0.05 1.5 5]);
plot(data_p(1, :),data_p(5,:),"b;func;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(2,:),"k;P_3;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(3,:),"g;P_5;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(4,:),"r;P_7;", 'LineWidth',1.25);
legend(legend("show"), "location", "northwest");

subplot(1,2,2);
hold on;
grid on;
title('Графики функции и приближений сплайнами');
axis([a-0.05,b+0.05 2.5 5]);
plot(data_s(1, :),data_s(5,:),"b;func;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(2,:),"k;S_3;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(3,:),"g;S_5;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(4,:),"r;S_7;", 'LineWidth',1.25);
legend(legend("show"), "location", "northwest");

figure;
subplot(1,2,1);
hold on;
grid on;
title('Графики ошибки приближений полиномами');
axis([a-0.05,b+0.05]);
plot(data_p(1, :),data_p(2,:)-data_p(5,:),"k;P_3;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(3,:)-data_p(5,:),"g;P_5;", 'LineWidth',1.25);
plot(data_p(1, :),data_p(4,:)-data_p(5,:),"r;P_7;", 'LineWidth',1.25);

subplot(1,2,2);
hold on;
grid on;
title('Графики ошибки приближений сплайнами');
axis([a-0.05,b+0.05]);
plot(data_s(1, :),data_s(2,:)-data_p(5,:),"k;S_3;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(3,:)-data_p(5,:),"g;S_5;", 'LineWidth',1.25);
plot(data_s(1, :),data_s(4,:)-data_p(5,:),"r;S_7;", 'LineWidth',1.25);
file = fopen([dir,'ngpoly.bin'],'rb');
data_p = fread(file,Inf,'double');
fclose(file);
file = fopen([dir,'ngspline.bin'],'rb');
data_s = fread(file,Inf,'double');
fclose(file);


figure;
hold on;
grid on;
title('Зависимость ошибки от числа точек для функции с разрывом производной');
semilogy(3:size(data_p)(1)+2,data_p,'k;poly;');
semilogy(3:size(data_s)(1)+2,data_s,'b;spline;');

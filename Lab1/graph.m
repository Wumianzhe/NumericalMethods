pkg load symbolic
file = fopen('~/Documents/Numerical/build/lab1/polysmooth.bin','rb');
data = fread(file,[5,Inf],'double');
fclose(file);

a = 1;
b = 2.5;
c = 1.7;

f = @(x) cot(x) + x^2;
g = @(x) f(c) + abs(f(x) - f(c));
l = @(x) 2*x - log10(x);

syms x;
ff = f(x);
gg = g(x);
% сетка
function ret=meshGen(a,b,n,f)
  for i = 0:n
    ret(1,i+1) = a + (b-a)/n*i;
    ret(2,i+1) = f(ret(1,i+1));
  endfor
endfunction

% корневой полином
function ret=omega(m,x)
  ret = ones(1,size(x)(2));
  for x_i = m(1,:)
    ret .*= (x-x_i);
  endfor
endfunction

figure;
hold on;
grid on;
title('Графики гладкой функции и приближений');
axis([a-0.05,b+0.05 1.5 5]);
plot(data(1, :),data(5,:),"b;func;", 'LineWidth',1.25);
plot(data(1, :),data(2,:),"k;P_3;", 'LineWidth',1.25);
plot(data(1, :),data(3,:),"g;P_5;", 'LineWidth',1.25);
plot(data(1, :),data(4,:),"r;P_7;", 'LineWidth',1.25);

mesh3 = meshGen(a,b,3,f);
plot(mesh3(1,:),mesh3(2,:),"ok");
mesh5 = meshGen(a,b,5,f);
plot(mesh5(1,:),mesh5(2,:),"og");
mesh7 = meshGen(a,b,7,f);
plot(mesh7(1,:),mesh7(2,:),"or");

figure;
hold on;
grid on;
title('Ошибка для приближений гладкой функции');
axis([a,b]);
plot(data(1, :),data(2,:)-data(5,:),"k;\\DeltaP_3;", 'LineWidth',1.25);
plot(data(1, :),data(3,:)-data(5,:),"g;\\DeltaP_5;", 'LineWidth',1.25);
plot(data(1, :),data(4,:)-data(5,:),"r;\\DeltaP_7;", 'LineWidth',1.25);

## ffd = diff(ff,x,8);
## df = function_handle(ffd);
## theor = @(x) abs(df(1))*abs(omega(mesh7,x))/factorial(8);
## plot(data(1, :),theor(data(1,:)), "b;theor;", 'LineWidth', 1.25);

file = fopen('~/Documents/Numerical/build/lab1/pointsmooth.bin','rb');
data = fread(file, Inf, 'double');
fclose(file);

figure;
hold on;
grid on;
title('Ошибка в зависимости от числа точек')
semilogy(3:size(data)(1)+2,data,'LineWidth',1.25);

file = fopen('~/Documents/Numerical/build/lab1/polycorner.bin','rb');
data = fread(file,[5,Inf],'double');
fclose(file);

figure;
hold on;
grid on;
title('Графики функции с углом и приближений');
axis([0.95,2.55 1.5 5]);
plot(data(1, :),data(5,:),"b;func;", 'LineWidth',1.25);
plot(data(1, :),data(2,:),"k;P_3;", 'LineWidth',1.25);
plot(data(1, :),data(3,:),"g;P_5;", 'LineWidth',1.25);
plot(data(1, :),data(4,:),"r;P_7;", 'LineWidth',1.25);

mesh3 = meshGen(a,b,3,g);
plot(mesh3(1,:),mesh3(2,:),"ok");
mesh5 = meshGen(a,b,5,g);
plot(mesh5(1,:),mesh5(2,:),"og");
mesh7 = meshGen(a,b,7,g);
plot(mesh7(1,:),mesh7(2,:),"or");

figure;
hold on;
grid on;
title('Ошибка для приближений функции с углом');
axis([1,2.5]);
plot(data(1, :),data(2,:)-data(5,:),"k;\\DeltaP_3;", 'LineWidth',1.25);
plot(data(1, :),data(3,:)-data(5,:),"g;\\DeltaP_5;", 'LineWidth',1.25);
plot(data(1, :),data(4,:)-data(5,:),"r;\\DeltaP_7;", 'LineWidth',1.25);

## ggd = diff(gg,x,8);
## dg = function_handle(ggd);
## theor = @(x) abs(dg(1))*abs(omega(mesh7,x))/factorial(8);
## plot(data(1, :),theor(data(1,:)), "c;theor;", 'LineWidth', 1.25);

file = fopen('~/Documents/Numerical/build/lab1/pointcorner.bin','rb');
data = fread(file, Inf, 'double');
fclose(file);

figure;
hold on;
grid on;
title('Ошибка для функции с углом в зависимости от числа точек');
semilogy(3:size(data)(1)+2,data,'LineWidth',1.25);

figure;
grid on;
title('Расстояние от точки разрыва до ближайшей точки сетки')
dist = @(n) min(abs(linspace(1,2.5,n+1) - 1.7));
for i = 3:50
  r(i-2) = dist(i);
endfor
plot(3:50,r);

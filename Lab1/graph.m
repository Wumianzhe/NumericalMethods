% sample size for each method
size = 14;
file = fopen("results.csv","r");
format = '%e, %f, %e, %d';
poly_bin = fscanf(file, format, [4,size]);
transc_bin = fscanf(file, format, [4,size]);
poly_sec = fscanf(file, format, [4,size]);
transc_sec = fscanf(file, format, [4,size]);
poly_fzero = fscanf(file, format, [4,size]);
transc_fzero = fscanf(file, format, [4,size]);
fclose(file);

% replacing all 0 with 1e-16
my_abs = @(x) max(abs(x),1e-16);

figure; axis([1,14]);
hold on;
grid on;
semilogy(1:14,my_abs(poly_bin(3,1:14)),'r;bin;','LineWidth',1.5);
semilogy(1:14,my_abs(poly_sec(3,1:14)),'g;secant;','LineWidth',1.5);
semilogy(1:14,my_abs(poly_fzero(3,1:14)),'b;fzero;','LineWidth',1.5);
semilogy(1:14,poly_bin(1,1:14),'k;eps;','LineWidth',1.25);
xlabel('eps');
ylabel('|f(x_e)|');
title('Polynom divergence');

figure;
axis([1,14]);
hold on;
grid on;
plot(1:14,poly_bin(4,1:14),'r;bin;','LineWidth',1.5);
plot(1:14,poly_sec(4,1:14),'g;secant;','LineWidth',1.5);
plot(1:14,poly_fzero(4,1:14),'b;fzero;','LineWidth',1.5);
xlabel('eps');
ylabel('n(e)');
title('Polynom iterations');

figure;
axis([1,14]);
hold on;
grid on;
semilogy(1:14,my_abs(transc_bin(3,1:14)),'r;bin;','LineWidth',1.5);
semilogy(1:14,my_abs(transc_sec(3,1:14)),'g;secant;','LineWidth',1.5);
semilogy(1:14,my_abs(transc_fzero(3,1:14)),'b;fzero;','LineWidth',1.5);
semilogy(1:14,poly_bin(1,1:14),'k;eps;','LineWidth',1.25);
xlabel('eps');
ylabel('|f(x_e)|');
title('Transcendental divergence');

figure;
axis([1,14]);
hold on;
grid on;
plot(1:14,transc_bin(4,1:14),'r;bin;','LineWidth',1.5);
plot(1:14,transc_sec(4,1:14),'g;secant;','LineWidth',1.5);
plot(1:14,transc_fzero(4,1:14),'b;fzero;','LineWidth',1.5);
xlabel('eps');
ylabel('n(e)');
title('Transcendental iterations');

figure;
title('Visual representation');
hold on;
axis([0.2,1.8]);
% using secant method to solve polynom with eps=1e-3
M_2 = 88; m_1=2;
cur=1.85; next=1.8;
my_eps=1e-3;
f = @(x) 2*x.^4+8*x.^3+8*x.^2-1;
x = 0.1:0.001:1.8;
plot(x,f(x));

while (abs(power(abs(next-cur),1.61)*M_2/2/m_1) > my_eps)
  prev=cur;
  cur=next;
  next=cur+(cur-prev)/(f(prev)/f(cur)-1);
  line([prev,next],[f(prev),0],'LineWidth',1.1);
  line([next,next],[f(next),0],'LineStyle','--','LineWidth',1.1);
end
line([0.2,1.8],[0,0],'Color','k','LineWidth',1.2);
plot(next,0,'r0','markersize',12)
grid on;

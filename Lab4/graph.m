file = fopen('~/Documents/Numerical/build/4-7.csv');
data_7 = fscanf(file, '%e, %e, %d, %d, %e, %e', [6,Inf]);
fclose(file);

figure;
hold on;
grid on;
loglog(data_7(1, :), data_7(2, :), "r;eigval;", 'LineWidth', 1.25);
loglog(data_7(1, :), data_7(5, :), "b;eigvec;", 'LineWidth', 1.25);
loglog(data_7(1, :), data_7(6, :), "g;discrepancy;", 'LineWidth', 1.25);
loglog(data_7(1, :), data_7(1, :), "k;eps;");

figure;
hold on;
grid on;
semilogx(data_7(1, :), data_7(3, :), "r;Jacobi;", 'LineWidth', 1.25);
semilogx(data_7(1, :), data_7(4, :), "b;INVIT;", 'LineWidth', 1.25);

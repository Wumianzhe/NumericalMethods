file = fopen("~/Documents/Numerical/build/3-7.csv");
data_7 = fscanf(file, '%e, %d, %e, %e', [ 4, inf ]);
fclose(file);

figure;
hold on;
grid on;
semilogx(data_7(1, :), data_7(2, :), "r;iterations;", 'LineWidth', 1.25);

figure;
hold on;
grid on;
loglog(data_7(1, :), data_7(3, :), "r;error;", 'LineWidth', 1.25);
loglog(data_7(1, :), data_7(4, :), "b;discrepancy;", 'LineWidth', 1.25);
loglog(data_7(1, :), data_7(1, :), "k;eps;")

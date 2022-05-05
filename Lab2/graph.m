file = fopen("~/Documents/Numerical/build/2-cond.csv");
data_cond = fscanf(file, '%e, %e, %e', [ 3, 30 ]);
fclose(file);
file = fopen("~/Documents/Numerical/build/2-dist.csv");
data_dist = fscanf(file, '%e, %e, %e', [ 3, 35 ]);
fclose(file);

figure;
hold on;
grid on;
loglog(data_cond(1, 1 : 30), data_cond(2, 1 : 30), 'r;ошибка;', 'LineWidth', 1.5)
loglog(data_cond(1, 1 : 30), data_cond(3, 1 : 30), 'g;невязка;', 'LineWidth', 1.5)

figure;
hold on;
grid on;
loglog(data_dist(1, :), data_dist(2, :), 'g;good;', 'LineWidth', 1.5);
loglog(data_dist(1, :), data_dist(3, :), 'r;bad;', 'LineWidth', 1.5);
xlabel('||\delta b||/||b||')

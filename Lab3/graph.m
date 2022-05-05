file = fopen("~/Documents/Numerical/build/lab3/results.csv");
data = fscanf(file, "%d, %e, %e", [3, Inf]);
fclose(file);

a = -2.0;
b = -0.3;

figure;
hold on;
grid on;
loglog(data(3,:), data(2,:), "k;error;",'LineWidth', 1.25);
loglog(data(3,:), data(3,:), "b;eps;");
xlabel("eps");

figure;
hold on;
grid on;
semilogx(data(3,:), data(1,:), "k;divisions;",'LineWidth', 1.25);
xlabel("eps");

file = fopen("~/Documents/Numerical/build/lab3/steps.csv");
data = fscanf(file, "%d, %e", [2, Inf]);
fclose(file);

figure;
hold on;
grid on;
h = (b-a)./(2.^data(1,:));
loglog(h, data(2,:), "k;error;", 'LineWidth', 1.25);
loglog(h,h.^2, "b;h^2;")
xlabel("h")

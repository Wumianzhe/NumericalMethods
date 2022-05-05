lambda = sort(eig(A));
alpha = 2/(lambda(1) + lambda(end));
C = eye(10) - alpha * A;

% clear;
f = fopen('results.csv','a');
def = optimset('fzero');
poly = @(x) 2*x^4 + 8*x^3 + 8*x^2 -1;
my_eps = 1e-1;
for i=1:14
  options=optimset(def,'TolX',my_eps);
  [res, fval, info, out]=fzero(poly, [1/9,1.840], options);
  fprintf(f,'%e, %0.16f, %0.12e, %d\n', my_eps, res, fval, out.iterations);
  my_eps/=10;
end

transc = @(x) (x-3)*cos(x)-1;
my_eps = 1e-1;
for i=1:14
  options=optimset(def,'TolX',my_eps);
  [res, fval, info, out]=fzero(transc, [5,5.4], options);
  fprintf(f,'%e, %0.16f, %0.12e, %d\n', my_eps, res, fval, out.iterations);
  my_eps/=10;
end
fclose(f);

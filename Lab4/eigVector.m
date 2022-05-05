F = fopen("a.bin","r");
A = fread(F,Inf,"double");
fclose(F);
side = sqrt(length(A));
A = reshape(A,[side,side]);
F = fopen("a.bin","w");
[V, d] = eigs(A,1,'sm');
fwrite(F,V,"double");
fclose(F);

iiMax = 10;
ratio = zeros(1,iiMax);
for ii = 1:iiMax
rng(ii)
n = 1e3;
b = 16;
A = randn(n,n);
x = randn(n,1);
X = randn(n,b);


ratio(ii) = timeit(@() A*X)/timeit(@() A*x);
end
ratio = sort(ratio);

mean(ratio(2:iiMax-1))
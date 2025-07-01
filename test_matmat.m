rng(1)
n = 1e5;
b = 8;
A = sprandsym(n,1e-3);
x = randn(n,1);
X = randn(n,b);


timeit(@() A*X)/timeit(@() A*x)
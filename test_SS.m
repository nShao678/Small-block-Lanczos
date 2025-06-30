rng(1)
n = 300;
NMax = 7;
bMax = 4;
err = zeros(NMax,bMax);
m = 2^(bMax-1);
a1 = linspace(-3,-2,n)';
a = linspace(-0.5,0.5,m)';
a2 = linspace(2,3,n)';
n = 2*n+m;
A = [a1;a;a2];
for iterb = 1:bMax

    b = 2^(iterb-1);
    d = m/b;
    V = randn(n,b);
    for iterN = 1:NMax
        N = 8*iterN;
        if N>=d
            lambda = blockSS(A, N, V,d);
            err(iterN,iterb) = max(min(abs(lambda-a')));
        end
    end
    
end
save('dataCIRR')
function lambda = blockSS(A, N, V,d)
% A is a column vector containing all eigenvalues
[n,b] = size(V);
theta = linspace(0, 2*pi, N+1); theta(end) = [];
z = exp(1i*theta);
% plot(real(z),imag(z),'x')
w = z/N;
Q = [];

for k = 0:d-1
    Q0 = zeros(n,b);
    for j = 1:N
        R = V./(z(j) - A);
        Q0 = Q0 + w(j)*z(j)^k*R;
    end
%     coef = zeros(n,1);
%     for j = 1:N
%         coef = coef+w(j)*z(j)^k./(z(j) - A);
%     end
%     coef = real(coef);
%     norm(coef(301:364)-A(301:364).^k)
%     Q0 = V.*coef;
%     Q0 = orth(real(Q0));
%     norm(Q0,'fro')-norm(Q0(201:200+d,:),'fro')
%     if k>0
%         Q0 = Q0-Q*(Q'*Q0);
%         Q0 = Q0-Q*(Q'*Q0);
%     end
%     norm(Q0)
%     if norm(Q0)<1e-10
%         1
%     end
    Q = orth([Q,real(Q0)]);
end

% 
% H0 = zeros(b*d, b*d);
% H1 = zeros(b*d, b*d);
% 
% for i = 1:d
%     for j = 1:d
%         block0 = M(:,:,i+j-1);
%         block1 = M(:,:,i+j);
%         H0((i-1)*b+1:i*b,(j-1)*b+1:j*b) = block0;
%         H1((i-1)*b+1:i*b,(j-1)*b+1:j*b) = block1;
%     end
% end
T = Q'*(Q.*A);
T = (T+T')/2;
lambda = eig(T);
lambda = lambda(abs(lambda)<1);
lambda = sort(lambda,'ascend');



end

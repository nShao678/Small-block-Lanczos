function Q = Blanczos(A,B,m)

[n,nB] = size(B);
Q = zeros(n,m);
d = m/nB;
% alpha = cell(m,1);
% beta = cell(m,1);
[B,~] = qr(B,0);
% tol = 1e-4;
for j = 1:d
    Q(:,(j-1)*nB+1:j*nB) = B;
    Z = A(B);
    % alpha{j} = B'*Z;
    Z = Z-Q(:,1:j*nB)*(Q(:,1:j*nB)'*Z);
    Z = Z-Q(:,1:j*nB)*(Q(:,1:j*nB)'*Z);
    % [B,beta{j}] = qr(Z,0);
    [B,~] = qr(Z,0);
end
% alpha = alpha(1:m);
% beta = beta(1:m-1);
% T = blkdiag(alpha{:});
% for ii = 1:m-1
%     T((ii-1)*nB+1:ii*nB,ii*nB+1:(ii+1)*nB) = beta{ii}';
%     T(ii*nB+1:(ii+1)*nB,(ii-1)*nB+1:ii*nB) = beta{ii};
% end
% T = (T+T')/2;
% nz = sum(abs(eig(T))<tol); 
end

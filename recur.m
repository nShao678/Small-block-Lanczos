function [S,Omegah] = recur(A,B,Omega,ii)


[b,d] = size(A);
idx = [ii,1:ii-1,ii+1:d];
A = A(:,idx);
B = B(:,idx);
Omega = Omega(:,idx);

Omegah = cell(1,d);
Bh = cell(1,d);
for ii = d:-1:1
    S = eye(b);
    for jj = ii+1:d
        S = B{ii}*S-S*Bh{jj};
    end
    Omegah{ii} = Omega{ii}*S;
    Bh{ii} = Omegah{ii}\(Omegah{ii}.*A(:,ii));
end

end


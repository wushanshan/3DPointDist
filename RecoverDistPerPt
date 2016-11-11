function [Dout,P] = RecoverDistPerPt(D, noise)
%D is the given n-by-n Euclidean distance matrix, we will add Guassian noise to it
%Dout is the recovered distance matrix, P is a n(n-1)-by-n(n-1) 0/1 matirx
%that satisfies: 1. each row/coln sum = 1; 2. every entry>=0; 3. diagonal
%block is all zeros; 4. off-diagonal (n-1)-by-(n-1) blocks sum to 1
%--------add noise to D, and transform it to n(n-1)-by-n(n-1) matrix------------%
n = size(D,1);
m = n*(n-1);
D = D + sqrt(noise)*randn(n,n);
for i = 1:n
    D(i,i)=0;
end
Dt = reshape(D, 1, []);
nonzeros = find(Dt); % indices of nonzeros
Dt = Dt(nonzeros); % remove the zeros
Dt = abs(diag(Dt)*ones(m,m) - ones(m,m)*diag(Dt));
%--------optimize over P matrix------------%
T = []; %n(n-1)-by-n matrix, coln-i has n-1 consecutive 1s starting from row-((i-1)*(n-1)+1)
for i = 1:n
    t = zeros(1,n);
    t(i) = 1;
    T = [T; repmat(t,n-1,1)];
end
cvx_begin
    variable P(m,m) symmetric;
    minimize trace(Dt*P);
    subject to
        P*ones(m,1) == ones(m,1);
        ones(1,m)*P == ones(1,m);
        P >= 0;
        T'*P*T == ones(n,n)-eye(n);
cvx_end
%--------recover the distance matrix------------%
Dout = D;
Perm = T'*P;
for i = 1:n
    t1 = Dout(i,:)';
    temp = Perm(:,((i-1)*(n-1)+1):i*(n-1))*t1(find(t1));
    Dout(i,:) = temp';
end

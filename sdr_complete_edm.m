function [EDM, X] = sdr_complete_edm(D, lambda)
% Given a noisy distance matrix D, find the best fitting point
% configuration using semidefinite relaxation.
% Code adapted from https://arxiv.org/pdf/1502.07541v2.pdf
n = size(D, 1);
x = -1/(n + sqrt(n));
y = -1/sqrt(n);
V = [y*ones(1, n-1); x*ones(n-1) + eye(n-1)];
e = ones(n, 1);

cvx_begin sdp
 variable G(n-1, n-1) symmetric;
 B = V*G*V';
 E = diag(B)*e' + e*diag(B)' - 2*B;
 maximize trace(G)-lambda * norm(E - D, 'fro');
 subject to
   G >= 0;
cvx_end

[U, S, V] = svd(B);
EDM = diag(B)*e' + e*diag(B)' - 2*B;
X = sqrt(S)*V';

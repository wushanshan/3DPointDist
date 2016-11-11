function main_RecNoise(n)
%---------------generate n random points in d-dimension-------------%
d = 3;
X = randn(n,d);
G = X*X';
 % Euclidean distance matrix, note that we use sqrt, so its rank is not 5
D0 = sqrt(ones(n,1)*diag(G)' + diag(G)*ones(1,n) - 2*G);
[~,P0] = RecoverDistPerPt(D0, 0); % P0 is the groundtruth
%---------------find the minimum distance------------------%
Dt = reshape(D0, 1, []);
nonzeros = find(Dt); % indices of nonzeros
Dt = Dt(nonzeros); % remove the zeros
dmin = min(Dt);
dmax = max(Dt);
%----------------adding noise to the distance matrix--------------%
Pdist = []; % distance between recovered P matrix and the groundtruth
Ddist = []; % distance between recovered distance matrix and the groudtruth
DeNdist = []; % distance between denoised distance matrix and the groundtruth
J = eye(n) - ones(1,1)/n; 
noise = 0:0.1:ceil(dmax/dmin);
for i = 1:length(noise);
    [Dout,P] = RecoverDistPerPt(D0, noise(i)*dmin); 
    Pdist = [Pdist, norm(P-P0,'fro')];
    Ddist = [Ddist, norm(Dout-D0,'fro')/norm(D0,'fro')];
    [DNout,~] = sdr_complete_edm(Dout.^2, 1); % TODO: need to tune the lambda parameter (currently set it as 1)
    DeNdist = [DeNdist, norm(D0-sqrt(DNout),'fro')/norm(D0,'fro')];
end
figure;
plot(noise, Pdist);
xlabel('noise/dmin')
ylabel('||Pout-P||_F')
figure;
plot(noise, Ddist);
hold on
plot(noise, DeNdist);
xlabel('noise/dmin')
ylabel('||Dout-D||_F/||D||_F')
legend('Recovered', 'Recovered and denoised')

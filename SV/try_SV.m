clear all
close all

% theta = [1, 0.97, 0.15^2];
beta = 0.05;
mu = 2*log(beta);
theta = [mu, 0.98, 0.2^2];

T = 1000;
[y,h_true] = generate_SV(theta,T);
plot(y)
plot(h_true)

mu = theta(1);
phi = theta(2);
sigma2 = theta(3);

% P0 = sqrt(sigma2/(1-phi^2)); % unconditional st. dev. 0.6120
% bin_range = 5*P0;
bin_range = 4;
N_bin = 30; %50;
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities
stdev_y = exp((mu+bin_midpoint)/2);

delta.h = 0.1;

h_init = var(y)*ones(1,T); 
h = h_init;

M = 10000;
H = zeros(M,T);
accept = zeros(M,1);
A_sum = zeros(M,1);

tic
% figure(1290)
% plot(h_true,'r')
% hold on
for ii = 1:M
    [h, acc, A_s] = update_h(y-mean(y),h, theta, delta.h);
    if (mod(ii,1000) == 1)
        toc;
%         plot(h,'color',[0 (ii/M) (1-ii/M)])
    end
    H(ii,:) = h;
    accept(ii,1) = acc;
    A_sum(ii,1) = A_s;
end
% hold off
toc

figure(120)
hold on
plot(h)
plot(h_true,'r')
plot(mean(H((5000+1):M,:),1),'g')
hold off
prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];

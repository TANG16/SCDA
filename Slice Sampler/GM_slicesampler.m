clear all

A = 1; B = 0; C1 = 3; C2 = 3; 
L = true;
mu_init = [0, 0.1];
kernel = @(a) GelmanMeng(a,A,B,C1,C2,L);


N = 2000;
Sample = slicesample(mu_init,N,'logpdf',kernel,'thin',1,'burnin',1000);  
%  [rnd, neval] = slicesample(initial=mu_init;nsamples=N;varargin={'logpdf',kernel,'thin',1,'burnin',1000})


% Plot    
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
 
x = -1:0.05:6;
n = length(x);
[X1,X2] = meshgrid(x,x);
V1 = reshape(X1,n*n,1); V2 = reshape(X2,n*n,1);
V = [V1,V2];

GM_fix = @(a) GelmanMeng(a,A,B,C1,C2,false);
GM = arrayfun(@(ii) GM_fix(V(ii,:)), 1:n*n, 'un', 0);
GM = reshape(GM,n,n);
GM = cell2mat(GM);
subplot(1,2,1)
GM_surf = surf(x,x,GM);
set(GM_surf,'LineStyle','none')  
subplot(1,2,2)
scatter(Sample(:,1),Sample(:,2));


%% ARCH
T = 1000;

sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
sigma1_k = sigma1/sqrt(kappa);
sigma2_k = sigma2/sqrt(kappa);

mu2 = 0; % gama = mu - mu2;
omega = 1;
alpha = 0.1;

eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
eps = eps/sqrt(kappa);  
y = zeros(T,1);
h_true = zeros(T,1);
h_true(1,1) = omega;
y(1,1) = sqrt(h_true(1,1)).*eps(1,1);
for ii = 2:T
    h_true(ii,1) = omega*(1-alpha) + alpha*(y(ii-1,1)).^2;
    y(ii,1) = sqrt(h_true(ii,1)).*eps(ii,1);
end
 y_S = var(y(1:T));
  
kernel_init = @(xx) -posterior_arch1_mex(xx, y(1:T), y_S)/T;
kernel = @(xx) posterior_arch1_mex(xx, y(1:T), y_S);

mu_init = [0, 1, 0.1, 0.05];
Sample = slicesample(mu_init,N,'logpdf',kernel,'thin',1,'burnin',1000);  
ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
for ii = 1:4
    subplot(2,2,ii)
    hist(Sample(:,ii),20)
end


ff = figure(3);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
for ii = 1:4
    subplot(2,2,ii)
    autocorr(Sample(:,ii),40)
end


ff = figure(4);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
for ii = 1:4
    subplot(2,2,ii)
    plot(Sample(:,ii))
end



%% Xi'an Blog: Beta(0.5,0.5)
kernel = @(xx) betapdf(xx,0.5,0.5);
figure(6)
plot(0:0.01:1,betapdf(0:0.01:1,0.5,0.5))
hold on
plot(0:0.01:1,betapdf(0:0.01:1,1,1),'m')
hold off

init = rand;
Sample = slicesample(init,N,'pdf',kernel,'thin',1,'burnin',1000);  
[hh,xs] = hist(Sample,50);
hh = hh/length(Sample)
figure(6)
hold on
bar(xs,hh,1)
set(gca,'XLim',[0 numel(hh)+1]);
plot(xs,betapdf(xs,0.5,0.5),'r')


figure(7)
hold on
plot(0:0.01:1,betapdf(0:0.01:1,0.5,0.5),'r')


dote = @(x,y) scatter(x,y,'MarkerFaceColor', [0.5 0 0.5]);
mote = @(x,y,z,w) line([x,z],[y,w]);

cst = betapdf(0.5,0.5,0.5)*0.5; %normalising constant
%inverting f(x)=d, 2nd degree equation
hitden= @(d) 0.5+0.5*sqrt(1-4*( cst/ max(d,betapdf(0.5,0.5,0.5)))^2).*[-1,1];
%output
% curve(dbeta(x,.5,.5),0,1,ylab='density',lwd=2,col='steelblue',n=1001);
x = rand;
u = rand*betapdf(x,0.5,0.5);
dote(x,u);
for t = 1:100 %100 slice steps
   bo = hitden(u);
%    nx=sample(c(runif(1,0,bo[1]),runif(1,bo[2],1)),1);
   el1 = bo(1)*rand;
   el2 = (1-bo(2))*rand + bo(2);
   if rand<0.5
       nx = el1;
   else
       nx = el2;
   end
   
   nu = rand*betapdf(nx,0.5,0.5);
   mote(x,u,nx,nu)
   x = nx;
   u = nu;
   dote(x,u)
end



%% try multislice
% https://nl.mathworks.com/matlabcentral/answers/276220-slice-sampling-for-multivariate-distribution

theta_init = [ones(1,12),0.5*ones(1,3)];
kernel = @(xx) try_multi_slice(xx);
W = [100*ones(1,12),0.5*ones(1,3)];
Sample = slicesample(theta_init,N,'logpdf',kernel,'width',W,'thin',1,'burnin',1000);
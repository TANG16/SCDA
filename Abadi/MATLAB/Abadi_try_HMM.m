clear all
close all

Abadi_Data_HMM

% v = randn(1,7);
% bp = randn(1,(ti-1));
% fec = 5*rand(1,ti-1);   %runif((ti-1),0,5),
v= [-2.7, 3.4, 0.2, 0.3, -1.6, -2.0, 1.3];
bp = 1.5*sin(1:ti-1);
fec = abs(cos(1:ti-1));

% Capture recapture data
[loglik_m, phij, phia, phijM, phiaM, lambda] = Abadi_Bugs_capture(v, mfem, mmal, bp, stdT);
loglik_Abadi_capture  = @(xx) -Abadi_Bugs_capture(xx, mfem, mmal, bp, stdT)/(2*numel(mfem));
theta_init = randn(1,5);
 
[theta_mle_cr, ~, ~, ~, ~, Sigma_mle_cr] = fminunc(loglik_Abadi_capture, theta_init);
Sigma_mle_cr = inv(2*numel(mfem)*Sigma_mle_cr);   
std_mle = sqrt(diag(Sigma_mle_cr));
% theta_mle_cr = [-2.6526    2.8457    0.3981    0.4366   -1.4739];


NadSurv = round(5+5*rand(1,ti)); %[NaN, round(5+5*rand(1,ti-1))];    %c(NA,round(runif(ti-1,5,10),0));
NadImm = round(5+5*rand(1,ti)); %[NaN, round(5+5*rand(1,ti-1))];  %c(NA,round(runif(ti-1,5,10),0));
% NadSurvprior = 20;
% Nadimmprior = 5;

Up = 100;
Na_prior = ones(1,Up+1)/(Up+1) ;
NadSurvprior = 20;
NadImmprior = 5;
% dummy = zeros(1,ti);
N_max= 20-1;
% C = 1000000;

%% Define the regression equations: best model structure (phi(a2+sex+T),p(sex+t),b(t)) for little owl
% Juvenile survival rate
index = v(1) + v(3) + v(4)*stdT;  % Male
phijM = exp(index)./(1+exp(index));  
index = v(1) + v(4)*stdT;          % Female
phij = exp(index)./(1+exp(index)); 
% Adult survival
index = v(1) + v(2) + v(3) +  v(4)*stdT;   % Male
phiaM = exp(index)./(1+exp(index));
index = v(1) + v(2) + v(4)*stdT;           % Female
phia =  exp(index)./(1+exp(index));
% Recapture rate
index = v(5) + bp;     % Male
lambdaM  = exp(index)./(1+exp(index));
index  = bp;	      % Female
lambda  = exp(index)./(1+exp(index));
% Immigration
index = v(6) + v(7)*voleH;     % Immigration rate as a function of 
im  = exp(index);              % vole abundance
                               % voleH is categorical with 2 levels (high,low)
  
%*******************************************
% Define the priors for the parameters %%%%
%*******************************************

pdf_v = normpdf(v, 0, sqrt(1/0.01));%I(-10,10)

% truncated
pdf_bp = normpdf(bp.*double((bp>-10) & (bp<10)), 0, sqrt(1/0.0001))/(normcdf(10,0, sqrt(1/0.0001))-normcdf(-10,0, sqrt(1/0.0001))); %I(-10,10)
prf_fec =  0.1 + 0*fec;% dunif(0,10)

%% %%%% The Integrated population model %%

%% Likelihood for reproductive data
rho = sample_size(1,1:ti-1).*fec;
pdf_nestlings = exp(-rho + nestlings(1,1:ti-1).*rho - sum(1:1:log(nestlings)));  %~ dpois(rho(ii))
  
%% Likelihood for population survey data 

% Due to the zeros trick: use a discrete uniform prior so the only influence on the posterior distr is the upper limit
% for t = 1:ti 
%     NadSurv(t) ~ dcat(Nad_prior[]); % Nad_prior = rep(1/(Up+1), Up+1); entered as data
%     NadImm(t) ~ dcat(Nad_prior[]) % Nad_prior = rep(1/(Up+1), Up+1); entered as data
% end
  
Ntot = zeros(N_max+1,ti);
G = NaN*zeros(N_max+1,N_max+1,ti);
P = zeros(N_max+1,ti);
Q = zeros(N_max+1,ti);
loglam = zeros(N_max+1,ti-1);
logmu = zeros(N_max+1,ti);


loglik = zeros(1,ti);

logfact = @(xx) sum(log(1:1:xx));

for t = 1:ti
    for jj = 0:N_max
        Ntot(jj+1,t) = jj + NadSurv(t) + NadImm(t); 
        Q(jj+1,t) = exp(-Ntot(jj+1,t) + popcount(t)*log(Ntot(jj+1,t)) - logfact(popcount(t)));
        logmu(jj+1,t) = log(im(t)) + log(Ntot(jj+1,t));
        if (t < ti)
            loglam(jj+1,t) = log(0.5) + log(fec(t)) + log(phij(t)) + log(Ntot(jj+1,t)); 
        end
    end
end

for jj = 0:N_max
    for ii = 0:N_max
      G(jj+1,ii+1,1) = 1/(N_max+1);
      G(jj+1,ii+1,2) = 1/(N_max+1);
    end 
    P(jj+1,1) = Nad_prior(1); % prior for Na(1)  observation probabilities
    
    if ((Ntot(jj+1,1) - NadSurv(2))>0)
        P(jj+1,2) = exp(NadSurv(2)*log(phia(1)) + (Ntot(jj+1,1) - NadSurv(2))*log(1-phia(1))... % part for survivors
        + logfact(Ntot(jj+1,1)) - logfact(Ntot(jj+1,1) - NadSurv(2)) - logfact(NadSurv(2))...
        - exp(logmu(jj+1,1)) + NadImm(2)*logmu(jj+1,1) - logfact(NadImm(2))); % part for immigrants
    else
        P(jj+1,2) = 0;
    end
end

 NadSurv(1) = NadSurvprior-1; % = 20
 NadImm(1) = NadImmprior-1; % = 5

 
for t = 3:ti
    % from 0 !!!
    for jj =  0:N_max  % old state N1_{t-1end= j
        for ii = 0:(N_max-1)  % new state N1_t=i
            % logfact is the log of the factorial: log(x!)
            G(jj+1,ii+1,t) = exp(-exp(loglam(jj+1,t-2)) + ii*loglam(jj+1,t-2) - logfact(ii));       
        end
        G(jj+1,(N_max+1),t) = max(0,1- sum(G(jj+1,1:(N_max),t)));
% G(i,j) = P(N1_{t-1}=j|N1_{t-2}=i)
% P(i,j) = P(Na_{t}|N1_{t-1}=i)
        if ((Ntot(jj+1,t-1) - NadSurv(t))>0)
            P(jj+1,t) = exp(NadSurv(t)*log(phia(t-1)) + (Ntot(jj+1,t-1) - NadSurv(t))*log(1-phia(t-1))... % part for survivors
            + logfact(Ntot(jj+1,t-1)) - logfact(Ntot(jj+1,t-1) - NadSurv(t)) - logfact(NadSurv(t))...
            - exp(logmu(jj+1,t-1)) + NadImm(t)*logmu(jj+1,t-1) - logfact(NadImm(t))); % part for immigrants
        else
            P(jj+1,t) = 0;
        end
    end  
    loglik(t) = sum(sum(G(:,:,t) * (P(:,t) .* Q(:,t)))); 
end
   
figure(10)
plot(loglik)



figure(1)
hold all
for ii = 1:ti
    plot(Q(:,ii))
end
hold off

figure(2)
hold all
for ii = 1:ti
    plot(P(:,ii))
end
hold off
 


figure(4)
hold all
for ii = 1:N_max
    plot(G(ii,:,26))
end
hold off
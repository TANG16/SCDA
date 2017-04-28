function [loglik_KF, N_smooth, phi1, phia, rho] = BKM_Bugs_statespace(theta, y, f, stdT)
% loglik_y
    n = max(size(y));
    alpha1 = theta(:,1);
    alphaa = theta(:,2);
    alphar = theta(:,3);
    beta1 = theta(:,4);
    betaa = theta(:,5);
    betar = theta(:,6);
    sigy = theta(:,7);
 
    index = alpha1 + beta1*f;
    phi1 = exp(index)./(1+exp(index));  
    index = alphaa + betaa*f; 
    phia = exp(index)./(1+exp(index));  	      
    index = alphar + betar*stdT;
    rho = exp(index);
    
    KFS_params = BKM_SSM(phi1,phia,rho,sigy);
    [N_smooth, ~, nu, F_inv]  = BKM_KFS_multi(y, KFS_params);
    loglik_KF = - 0.5*(n*log(2*pi) - sum(log(abs(F_inv))) + sum(nu.*F_inv.*nu)); 
% 	% Define the initial population priors
% %     N1 = zeros(T,1);
% %     Na = zeros(T,1);
% 	
%     loglik_N1 = 0;
%     loglik_Na = 0;   
%     sig1 = 1000^2;
%     siga = 1000^2;
%     for t = 1:2
%         loglik_N1 = loglik_N1  - 0.5*(log(2*pi) + log(sig1) + ((N1(t)-200)^2)/(sig1));
%         loglik_Na = loglik_Na  - 0.5*(log(2*pi) + log(siga) + ((Na(t)-1000)^2)/(siga));
% 
% % 		N1(t)  = 200 +  1000*randn; %dnorm(200,0.000001)
% % 		Na(t)  = 1000 + 1000*randn; %dnorm(1000,0.000001)
%     end
% 
% 	% Define the system process for the census/index data using the Normal approximation
%     for t = 3:T
% 		mean1 = rho(t-1)*phi1(t-1)*Na(t-1);
% 		meana = phia(t-1)*(N1(t-1)+Na(t-1));
% 		
% % 		tau1(t) = 1/(Na(t-1)*rho(t-1)*phi1(t-1));
% % 		taua(t) = 1/((N1(t-1)+Na(t-1))*phia(t-1)*(1-phia(t-1)));
% 		sig1 = sqrt(Na(t-1)*rho(t-1)*phi1(t-1));
% 		siga = sqrt((N1(t-1)+Na(t-1))*phia(t-1)*(1-phia(t-1)));
% 		
% % 		N1(t) = mean1 + sig1*randn;
% % 		Na(t) = meana + siga*randn;
%         loglik_N1 = loglik_N1  - 0.5*(log(2*pi) + log(sig1) + ((N1(t)-mean1)^2)/(sig1));
%         loglik_Na = loglik_Na  - 0.5*(log(2*pi) + log(siga) + ((Na(t)-meana)^2)/(siga));
%     end
% 	
% 	% Define the observation process for the census/index data
% 	loglik_y = 0;
%     for t = 3:T
% 	    %y(t) ~ dnorm(Na(t),tauy)
%         loglik_y = loglik_y - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na(t))^2)/(sigy));
%     end
end
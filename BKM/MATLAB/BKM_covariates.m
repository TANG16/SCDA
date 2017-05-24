function [phi1, phia, rho, lambda] = BKM_covariates(theta,f,stdT)    
    alpha1 = theta(:,1);
    alphaa = theta(:,2);
    alphar = theta(:,3);
    alphal = theta(:,4);
    
    beta1 = theta(:,5);
    betaa = theta(:,6);
    betar = theta(:,7);
    betal = theta(:,8);
%     sigy = theta(:,9);
 
    index = alpha1 + beta1*f;
    phi1 = exp(index)./(1+exp(index));  
    index = alphaa + betaa*f; 
    phia = exp(index)./(1+exp(index));  	      
    index = alphar + betar*stdT;
    rho = exp(index);   
    index = alphal + betal*stdT;
    lambda = exp(index)./(1+exp(index));  
end
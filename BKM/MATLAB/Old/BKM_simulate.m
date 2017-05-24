function [y, N1, Na] = BKM_simulate(theta,f,mode)

    T = length(f);
    time = 1:T;                
    stdT = (time - mean(time))/std(time);
    
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
    

    N1 = zeros(1,T);
    Na = zeros(1,T);
    
    if ((nargin == 2) || strcmp(mode,'approx'))
        N1(1,1) = 200 + 10^2*randn;
        Na(1,1) = 1000 + 10^2*randn;

        for ii = 2:T
            mu1 = Na(1,ii-1)*rho(1,ii-1)*phi1(1,ii-1);
            N1(1,ii) = mu1 + sqrt(mu1)*randn;

            mua = (N1(1,ii-1)+Na(1,ii-1))*phia(1,ii-1);
            siga = mua*(1-phia(1,ii-1));
            Na(1,ii) = mua + sqrt(siga)*randn;            
        end
    elseif strcmp(mode,'exact')
        N1(1,1) = poissrnd(200); % ?
        Na(1,1) = binornd(2000,0.5); % ?
        
        for ii = 2:T
            lambda = Na(1,ii-1)*rho(1,ii-1)*phi1(1,ii-1);
            N1(1,ii) =  poissrnd(lambda);
            
            mua = (N1(1,ii-1)+Na(1,ii-1));
            Na(1,ii) = binornd(mua,phia(1,ii));      
        end        
    else
        error('Incorrect mode!');
    end
    y = Na + sqrt(sigy)*randn(1,T);
    
end

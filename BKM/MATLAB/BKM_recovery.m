function loglik_m = BKM_recovery(theta, m, f, stdT)
    [T1 , T2] = size(m);
    T2 = T2-1;
    
   
    alpha1 = theta(:,1);
    alphaa = theta(:,2);
%     alphar = theta(:,3);
    alphal = theta(:,4);
    
    beta1 = theta(:,5);
    betaa = theta(:,6);
%     betar = theta(:,7);
    betal = theta(:,8);
    sigy = theta(:,9);
    
	% Define the logistic regression equations   
    index = alpha1 + beta1*f;
    phi1 = exp(index)./(1+exp(index));  
    index = alphaa + betaa*f; 
    phia = exp(index)./(1+exp(index));  	      
    index = alphal + betal*stdT;
    lambda = exp(index)./(1+exp(index));  
	 
	% Calculate the cell probabilities for the recovery table 
    % 	p = zeros(T1,T1);
    % Calculate the diagonal (dies and gets recovered immediately)
    p = diag(lambda.*(1-phi1));	      

    % Calculate value one above the diagonal
    indx = sub2ind(size(p),1:T1-1,2:T1);
    p(indx) = lambda(2:T1).* phi1(1:T1-1).*(1-phia(2:T1));
    
    % Calculate remaining terms above diagonal
    for t1 = 1:(T1-2)
        for t2 = (t1+2):T2
%             for t = (t1+1):(t2-1)
%                 lphi(t1, t2, t) = log(phia(t));
%             end
            lphi = sum(log(phia(t1:t2-2)));
            p(t1,t2) = lambda(t2)*phi1(t1)*(1-phia(t2))*exp(lphi);
%             p(t1,t2) = lambda(t2)*phi1(t1)*(1-phia(t2))*exp(sum(lphi(t1,t2,(t1+1):(t2-1))))
        end
    end

  	% Probability of an animal never being seen again
	p = [p, 1 - sum(p,2)];  % p(t1, T2+1) = 1 - sum(p(t1,1:T2));
    
%     logp = log(p);
%     logp(~isfinite(logp)) = 0;
%     loglik_m = sum(sum(m.*logp));
    loglik_m = sum(sum(log(mnpdf(m,p))));

%     loglik_m = 0;
%     for t = 1:T1
%         loglik_m = loglik_m + log(mnpdf(m(t,:),p(t,:)));
%     end
end
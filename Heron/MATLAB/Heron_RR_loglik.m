function loglik = Heron_RR_loglik(theta, m, f)
% function loglik = Heron_RR_loglik(m, lambda, phi1, phi2, phi3, phi4, T_ring)
% The ring-recovery data for individuals released between 1955 and 1997,
% i.e. between time indices t=t1,...,t2, where t1 = 28 and t2 = 70.
    T = length(f);
    T_ring = size(m,1);
    time = 1:T_ring;
    time = (time-mean(time))/std(time);
    year = 27 + (1:T_ring);
    
    alpha = theta(1:4);
    beta = theta(5:8);
    alphal = theta(9);
    betal = theta(10);

    ind = alpha(1) + beta(1)*f(1:T-1);
    phi1 = exp(ind)./(1+exp(ind));
    ind = alpha(2) + beta(2)*f(1:T-1);
    phi2 = exp(ind)./(1+exp(ind));
    ind = alpha(3) + beta(3)*f(1:T-1);
    phi3 = exp(ind)./(1+exp(ind));
    ind = alpha(4) + beta(4)*f(1:T-1);
    phi4 = exp(ind)./(1+exp(ind));

    ind = alphal + betal*(time(1:T_ring));
    lambda =  exp(ind)./(1+exp(ind));

    phi1 = phi1(year);
    phi2 = phi2(year);
    phi3 = phi3(year);
    phi4 = phi4(year);
    
    % Calculate the cell probabilities for the recovery table 
    % Calculate the diagonal
    p = diag(lambda.*phi1);	      

%     for t = 1:T_ring
%        p(t, t) = lambda(t) * (1-phi1(t))   ;
%     end

    % Calculate value one above the diagonal
    indx = sub2ind(size(p),1:T_ring-1,2:T_ring);
    p(indx) = lambda(2:T_ring).* phi1(1:T_ring-1).*(1-phi2(2:T_ring));
%     for t = 1:(T_ring-1)
%         p(t, t+1) = lambda(t+1) * phi1(t)*(1-phi2(t+1));
%     end

    % Calculate value two above the diagonal    
    indx = sub2ind(size(p),1:T_ring-2,3:T_ring);
    p(indx) = lambda(3:T_ring).* phi1(1:T_ring-2).*phi2(2:T_ring-1).*(1-phi3(3:T_ring));
%     for t = 1:(T_ring-2) 
%         p(t,t+2) = lambda(t+2) * phi1(t)*phi2(t+1)*(1-phi3(t+2));
%     end

  % Calculate remaining terms above diagonal 
    indx = sub2ind(size(p),1:T_ring-3,4:T_ring);
    p(indx) = lambda(4:T_ring).*phi1(1:T_ring-3).*phi2(2:T_ring-2).*phi3(3:T_ring-1).*(1-phi4(4:T_ring));
    for t1 = 1:(T_ring-4)
        for t2 = (t1+4):T_ring  
            lphi = sum(log(phi4((t1+3):(t2-1))));
            p(t1,t2) = lambda(t2)*phi1(t1)*phi2(t1+1)*phi3(t1+2)*(1-phi4(t2))*exp(lphi);
        end
    end
%     for t1 = 1:(T_ring-3)
%         for t2 = (t1+3):T_ring
%             for t = (t1+1):(t2-1)
%                 lphi(t1, t2, t) = log(phi4(t));
%             end
%             % Probabilities in table
%             p(t1,t2) = lambda(t2)*phi1(year(t1))*(1-phi4(year(t2)))*...
%                 exp(sum(lphi(t1,t2,(t1+1):(t2-1))));
%         end
%     end
% 
%     for(t1 = 1:T_ring) 
%         for(t2 = 1:(t1-1))
%             % Zero probabilities in lower triangle of table 
%             p(t1, t2) = 0
%         end
%         % Probability of an animal never being seen again 
%         p(t1, (T_ring+1)) = 1 - sum(p(t1, 1:T_ring))
%     end
 
    % Probability of an animal never being seen again
    p = [p, 1 - sum(p,2)];   
  
    logp = log(p);
    logp(~isfinite(logp)) = 0;
    loglik = sum(sum(m.*logp));
    
%     for t = 1:T_ring
%         m(t,1:(T_ring+1)) ~ dmulti(p(t, ), rel(t)) 
%     end
end
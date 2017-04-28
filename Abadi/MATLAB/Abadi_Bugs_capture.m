function [loglik_m, phij, phia, phijM, phiaM, lambda] = Abadi_Bugs_capture(v, mfem, mmal, bp, stdT)
    T1 = size(mfem,1);
    ti = T1/2;
    	
    % Define the logistic regression equations       
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

     % m-array cell probabilities for juveniles and adults
    q = 1 - lambda;
    qM = 1 - lambdaM;

    % Main diagonals 
    % females (juveniles and adults)
    pr = [diag(phij.*lambda);
          diag(phia.*lambda)];
    % males (juveniles and adults)
    prM = [diag(phijM.*lambdaM);
          diag(phiaM.*lambdaM)];

    %  Above main diagonal
    for ii = 1:(ti-1)  
        for jj = (ii+1):(ti-1)    
            pr(ii,jj) = phij(ii)*prod(phia((ii+1):jj))*prod(q(ii:(jj-1)))*lambda(jj);
            pr(ii+ti-1,jj) = prod(phia(ii:jj))*prod(q(ii:(jj-1)))*lambda(jj);
            prM(ii,jj) = phijM(ii)*prod(phiaM((ii+1):jj))*prod(qM(ii:(jj-1)))*lambdaM(jj);
            prM(ii+ti-1,jj) = prod(phiaM((ii+1):jj))*prod(qM(ii:(jj-1)))*lambdaM(jj);
        end 
    end 

    % Last column
    pr = [pr, 1- sum(pr,2)];
    prM = [prM, 1- sum(prM,2)];

    % for ii = 1:(2*(ti-1))
    %     m(ii,1:ti) ~ dmulti(pr(ii,:),r(ii))
    %     mM(ii,1:ti) ~ dmulti(prM(ii,:),rM(ii))
    % end 
    logpr = log(pr);
    logpr(~isfinite(logpr)) = 0;
    loglik_mfem = sum(sum(mfem.*logpr)); 

    logprM = log(prM);
    logprM(~isfinite(logprM)) = 0;
    loglik_mmal = sum(sum(mmal.*logprM)); 

    loglik_m = loglik_mfem + loglik_mmal;

%     loglik_m = 0;
%     for t = 1:T1
%         loglik_m = loglik_m + log(mnpdf(m(t,:),p(t,:)));
%     end
end
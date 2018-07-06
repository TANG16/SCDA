% SVRW.m
function [h, S] = SVRW(ystar,h,omega2h,Vh)
    T = length(h);
        %% parameters for the Gaussian mixture
    pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
    mui = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 ...
        -1.08819] - 1.2704; 
    sig2i = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
    sigi = sqrt(sig2i);
        %% sample S from a 7-point distrete distribution
    temprand = rand(T,1);
    q = repmat(pi,T,1).*normpdf(repmat(ystar,1,7),repmat(h,1,7) ... 
        +repmat(mui,T,1),repmat(sigi,T,1));
    q = q./repmat(sum(q,2),1,7);
    S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;
        %% sample h
    H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    invOmegah = spdiags([1/Vh; 1/omega2h*ones(T-1,1)],0,T,T);
    d = mui(S)'; invSigystar = spdiags(1./sig2i(S)',0,T,T);
    Kh = H'*invOmegah*H + invSigystar;
    Ch = chol(Kh,'lower');
    hhat = Kh\(invSigystar*(ystar-d));
    h = hhat + Ch'\randn(T,1);
end
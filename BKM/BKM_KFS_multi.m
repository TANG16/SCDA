function [alpha_smooth, V, nu, F_inv]  = BKM_KFS_multi(y, KFS_params)
% Univariate multistate Kalman Filter and Smoother with time varying parameters T and Q 
% Notation from DK 2012
% y_t = Z_t alpha_t + epsilon_t, epsilon_t ~ N(0, H_t)
% alpha_t+1 = T_t alpha_t + R_t eta_t, eta_t ~ N(0, Q_t)
% alpha_1 ~ N(a1, P1)
    n = max(size(y)); 
    m = size(KFS_params.P1,1);
    
%% Set parameters and initial values
    a = zeros(m,n+1);       % filtered state
    P = zeros(m,m,n+1);     % filtered state variance
    
    a(:,1) = KFS_params.a1;
    P(:,:,1) = KFS_params.P1;
    c = KFS_params.c;
    H = KFS_params.H;
    Q = KFS_params.Q;
    d = KFS_params.d;
    T = KFS_params.T;
    R = KFS_params.R;
    R1 = KFS_params.R1;
    R2 = KFS_params.R2;
    Z = KFS_params.Z;
    
%% Output of Kalman filter
 
    nu = zeros(1,n);            % prediction error
    F = zeros(1,n);             % prediction variance
    F_inv = zeros(1,n);
    K = zeros(m,n);             % = P/F, Kalman gain (regression coefficient of alpha on nu)
    L = zeros(m,m,n);
    
    for ii = 1:n
        nu(1,ii) = y(1,ii) - c - Z*a(:,ii);
        
        P_tmp = squeeze(P(:,:,ii));
        F(1,ii) = Z*P_tmp*Z' + H;
        F_inv(1,ii) = 1/F(1,ii);
        K(:,ii) = T(:,:,ii)*P_tmp*Z'*F_inv(1,ii);
        L(:,:,ii) = T(:,:,ii) - K(:,ii)*Z;
        
        a(:,ii+1) = d + T(:,:,ii)*a(:,ii) + K(:,ii)*nu(1,ii);
        R_tmp = [R1 * a(:,ii), R2 *  a(:,ii)];
        Q_tmp = R_tmp*squeeze(Q(:,:,ii));
        P(:,:,ii+1) = T(:,:,ii)*P_tmp*T(:,:,ii)' + R*Q_tmp*R' - ...
            T(:,:,ii)*P_tmp*Z'*K(:,ii)';
    end
    
%% State smoothing 
    alpha_smooth = zeros(m,n);  % smoothed state
    V = zeros(m,m,n);             % smoothed state variance
    r = zeros(m,n);             % smoothing cumulant
    N = zeros(m,m,n);             % smoothing variance cumulant

    for ii = n:-1:2
        r(:,ii-1) = Z'*F_inv(1,ii)*nu(1,ii) + L(:,:,ii)'*r(:,ii);
        alpha_smooth(:,ii) = a(:,ii) + P(:,:,ii)*r(:,ii-1);
        N(:,:,ii-1) = Z'*F_inv(1,ii)*Z + L(:,:,ii)'*N(:,:,ii)*L(:,:,ii);
        V(:,:,ii) = P(:,:,ii) - P(:,:,ii)*N(:,:,ii-1)*P(:,:,ii);
    end
    r0 = Z'*F_inv(1,1)*nu(1,1) + L(:,:,1)'*r(:,1);
    alpha_smooth(:,1) = a(:,1) + P(:,:,1)*r0;
    N0 = Z'*F_inv(1,1)*Z + L(:,:,1)'*N(:,:,1)*L(:,:,1);
    V(:,:,1) = P(:,:,1) - P(:,:,1)*N0*P(:,:,1);
    
%     theta_smooth = c + Z.*alpha_smooth;
    
%% Distrurbance smoothing   
%     eps_smooth = zeros(1,n);       % observation error
%     C = zeros(1,n);             % observation error variance
%     r = zeros(m,1);             % weighted sum of state innovations nu after t
%     % u - smoothing error 
%     N = zeros(m,m);              % state smoothing error variance
%     
%     for ii = n:-1:1
%         D = F_inv(1,ii) + K(:,ii)'*N*K(:,ii); % scalar
%         C(1,ii) = H - H*D*H; % scalar
%         u = F_inv(1,ii)*nu(1,ii) - K(:,ii)'*r;     % scalar 
%         eps_smooth(1,ii) = H*u;
%         r = Z'*u + T(:,:,ii)'*r; %Z'*nu(1,ii)*F_inv(1,ii) + r*squeeze(L(:,:,ii));
%         %N = Z*Z'*F_inv(1,ii) + squeeze(L(:,:,ii))*N*squeeze(L(:,:,ii));
%         N = Z'*D*Z + T(:,:,ii)'*N*T(:,:,ii) - Z'*K(:,ii)'*N*T(:,:,ii) -...
%             T(:,:,ii)'*N*K(:,ii)*Z; 
%     end
%     
% %% Recovering the smoothed signal: alpha_smooth = y-epsilon;
%     alpha_smooth = y - eps_smooth;
%     V = C;
%     
%     if (m == 1)
% %         P = squeeze(P);
%         L = squeeze(L);
%     end
end
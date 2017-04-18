clear all
close all

BKM_Data

%% Recovery Data
loglik_BKM_recovery = @(xx) -BKM_Bugs_recovery(xx, m, f_r, stdT_r)/numel(m);
theta_init = [1, 2, -4, -2, 0.1, -0.3];
% theta = [alpha1, alphaa, alphal, beta1, betaa, betal]

[theta_mle_r, ~, ~, ~, ~, Sigma_mle_r] = fminunc(loglik_BKM_recovery, theta_init);
Sigma_mle_r = inv(numel(m)*Sigma_mle_r);
% theta_mle = [0.5297    1.5241   -4.5708   -0.1659   -0.2937   -0.3452]
std_mle = sqrt(diag(Sigma_mle_r));
% std_mle = [0.0686    0.0699    0.0352    0.0617    0.0428    0.0390]
[~, phi1_r, phia_r] = BKM_Bugs_recovery(theta_mle_r, m, f, stdT_r);



%% Index Data
theta_init = [1, 2, -2, -2, 0.1, -0.7, 1];
theta = theta_init;
f = f_ss;
stdT = stdT_ss;
% theta = [alpha1, alphaa, alphar, beta1, betaa, betar, sigy]; 
% stdT = stdT_ss;
% f = f_ss;

loglik_BKM_statespace = @(xx) -BKM_Bugs_statespace(xx, y, f_ss, stdT_ss)/T;
options = optimoptions('fminunc');%,'MaxFunEvals',2000,'MaxIter',1500);
[theta_mle_ss, ~, ~, ~, ~, ~] = fminunc(loglik_BKM_statespace, theta_init,options);
[theta_mle_ss, ~, ~, ~, ~, ~] = fminunc(loglik_BKM_statespace, theta_mle_ss);
[theta_mle_ss, ~, ~, ~, ~, Sigma_mle_ss] = fminunc(loglik_BKM_statespace, theta_mle_ss);

[~, N_smooth, phi1_ss, phia_ss, rho] = BKM_Bugs_statespace(theta_mle_ss, y, f_ss, stdT_ss);




%% Combined
theta_init = [1, 2, -2, -4, -2, 0.1, -0.7, -0.3, 1];
% theta = [alpha1, alphaa, alphar, alphals, beta1, betaa, betar, betal, sigy]; 
theta = theta_init;
loglik_BKM_combined = @(xx) -BKM_Bugs_combined(xx, y, f, m)/(length(y)+numel(m));
% [loglik_C, loglik_m, loglik_KF, N_smooth, PHI1, PHIA, RHO, LAMBDA] = BKM_Bugs_combined(theta, y, f, m)
[theta_mle_c, ~, ~, ~, ~, Sigma_mle_c] = fminunc(loglik_BKM_combined, theta_init);
% [0.542 1.558 -1.161 -4.581 -0.163 -0.243 -0.355 -0.362 29218.134]
Sigma_mle_c = inv((numel(m)+length(y))*Sigma_mle_c);
% theta_mle = [0.5297    1.5241   -4.5708   -0.1659   -0.2937   -0.3452]
std_mle_c = sqrt(diag(Sigma_mle_c)) 
% [0.068 0.071 0.093 0.035 0.061 0.039 0.045 0.040 8013.534]


alpha1     0.5258 0.072
alphaa     1.4841 0.079
alphal    -4.5917 0.035 
alphar    -1.1017 0.094 
beta1     -0.1889 0.054 
betaa     -0.2360 0.038 
betal     -0.3574 0.039 
betar     -0.2871 0.037 
sigy   28822.1313 8817

[loglik_C, loglik_m, loglik_KF, N_smooth, PHI1, PHIA, RHO, LAMBDA] = BKM_Bugs_combined(theta_mle_c, y, f, m);

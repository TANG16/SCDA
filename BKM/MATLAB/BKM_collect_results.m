results_path = ['Results/BurnIn_',num2str(BurnIn),'/'];

name = [results_path,'BKM_DA_NRW_U.mat'];
load(name)
NN_DA = squeeze(NN(2,:,:));
Theta_DA = Theta;
time_sampl_DA = time_sampl;
mean_accept_DA = mean(accept(:,37:end));

name = [results_path,'BKM_adapt_Nq',num2str(10),'.mat'];
load(name)
NN_HMM_adapt_10 = NN;
Theta_HMM_adapt_10  = Theta;
time_sampl_HMM_adapt_10 = time_sampl;
mean_accept_HMM_adapt_10 = mean(accept);

name = [results_path,'BKM_adapt_Nq',num2str(20),'.mat'];
load(name)
NN_HMM_adapt_20 = NN;
Theta_HMM_adapt_20  = Theta;
time_sampl_HMM_adapt_20 = time_sampl;
mean_accept_HMM_adapt_20 = mean(accept);


name = [results_path,'BKM_adapt_Nq',num2str(30),'.mat'];
load(name)
NN_HMM_adapt_30 = NN;
Theta_HMM_adapt_30  = Theta;
time_sampl_HMM_adapt_30 = time_sampl;
mean_accept_HMM_adapt_30 = mean(accept);


name = [results_path,'BKM_bin_Nbin',num2str(10),'.mat'];
load(name)
NN_HMM_bin_10 = NN;
Theta_HMM_bin_10  = Theta;
time_sampl_HMM_bin_10 = time_sampl;
mean_accept_HMM_bin_10 = mean(accept);

name = [results_path,'BKM_bin_Nbin',num2str(20),'.mat'];
load(name)
NN_HMM_bin_20 = NN;
Theta_HMM_bin_20  = Theta;
time_sampl_HMM_bin_20 = time_sampl;
mean_accept_HMM_bin_20 = mean(accept);


name = [results_path,'BKM_bin_Nbin',num2str(30),'.mat'];
load(name)
NN_HMM_bin_30 = NN;
Theta_HMM_bin_30  = Theta;
time_sampl_HMM_bin_30 = time_sampl;
mean_accept_HMM_bin_30 = mean(accept);


clearvars -except -regexp ...
    '^NN_' '^Theta_' '^time_sampl_' 'mean_accept_' 'BurnIn'
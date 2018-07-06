% ind_meth = ones(8,1);
% BurnIn = 10000;
% ext = '_Selected';
results_path = ['Results/BurnIn_',num2str(BurnIn),ext,'/'];

if ind_meth(1)
    name = [results_path,'BKM_DA_NRW_U.mat'];
    load(name)
    NN_DA = squeeze(NN(2,:,:));
    Theta_DA = Theta;
    time_sampl_DA = time_sampl;
    mean_accept_DA = mean(accept(:,37:end));
end

if ind_meth(2)
    name = [results_path,'BKM_adapt_Nq',num2str(10),'.mat'];
    load(name)
    NN_HMM_adapt_10 = NN;
    Theta_HMM_adapt_10  = Theta;
    time_sampl_HMM_adapt_10 = time_sampl;
    mean_accept_HMM_adapt_10 = mean(accept);
end

if ind_meth(3)
    name = [results_path,'BKM_adapt_Nq',num2str(20),'.mat'];
    load(name)
    NN_HMM_adapt_20 = NN;
    Theta_HMM_adapt_20  = Theta;
    time_sampl_HMM_adapt_20 = time_sampl;
    mean_accept_HMM_adapt_20 = mean(accept);
end

if ind_meth(4)
    name = [results_path,'BKM_adapt_Nq',num2str(30),'.mat'];
    load(name)
    NN_HMM_adapt_30 = NN;
    Theta_HMM_adapt_30  = Theta;
    time_sampl_HMM_adapt_30 = time_sampl;
    mean_accept_HMM_adapt_30 = mean(accept);
end

if ind_meth(5)
    name = [results_path,'BKM_bin_Nbin',num2str(10),'.mat'];
    load(name)
    NN_HMM_bin_10 = NN;
    Theta_HMM_bin_10  = Theta;
    time_sampl_HMM_bin_10 = time_sampl;
    mean_accept_HMM_bin_10 = mean(accept);
end

if ind_meth(6)
    name = [results_path,'BKM_bin_Nbin',num2str(20),'.mat'];
    load(name)
    NN_HMM_bin_20 = NN;
    Theta_HMM_bin_20  = Theta;
    time_sampl_HMM_bin_20 = time_sampl;
    mean_accept_HMM_bin_20 = mean(accept);
end

if ind_meth(7)
    name = [results_path,'BKM_bin_Nbin',num2str(30),'.mat'];
    load(name)
    NN_HMM_bin_30 = NN;
    Theta_HMM_bin_30  = Theta;
    time_sampl_HMM_bin_30 = time_sampl;
    mean_accept_HMM_bin_30 = mean(accept);
end

if ind_meth(8)
%     name = [results_path,'BKM_exact_Nmax690.mat'];
    name = [results_path,'BKM_exact_Nmax679.mat'];
    load(name)
    NN_HMM_exact = NN;
    Theta_HMM_exact  = Theta;
    time_sampl_HMM_exact = time_sampl;
    mean_accept_HMM_exact = mean(accept);
end

clearvars -except -regexp ...
    '^NN_' '^Theta_' '^time_sampl_' 'mean_accept_' 'BurnIn' 'ext' 'ind_meth'
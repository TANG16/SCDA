clear all
close all
BurnIn = 10000;

ext = '';
% ext = '_Selected';
BKM_collect_results;

figures_path = ['Figures/BurnIn_',num2str(BurnIn),'/'];
%     figures_path = ['Figures/BurnIn_',num2str(BurnIn),'/'];
 
params = {'$\alpha_1$', '$\alpha_a$', '$\alpha_{\rho}$', '$\alpha_{\lambda}$', ...
            '$\beta_1$', '$\beta_a$', '$\beta_{\rho}$', '$\beta_{\lambda}$',...
            '$\sigma^{2}_{y}$'};
        
        
method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30', 'Exact'};
K = length(method);

TH_all = who('-regexp', '^Theta');
NN_all = who('-regexp', '^NN');
accept_all = who('-regexp','mean_accept');
time_sampl_all = who('-regexp','time');


% time
Time = zeros(8,1);
for jj = 1:8
    Time(jj) = eval(char(time_sampl_all{jj}));
end
Time_b20 = Time(6);
Time(6) = Time(7);
Time(7) = Time_b20;

% ESS    
T = 36; 
ESS_N_sig = zeros(8,T);
D = 9; 
ESS_TH_sig = zeros(8,D);

for jj = 1:K
    ESS_N_sig(jj,:) = ESS(eval(char(NN_all{jj}))',0);
    ESS_TH_sig(jj,:) = ESS(eval(char(TH_all{jj})),0);
end

% Na
N_ESS = ESS_N_sig;

[~, ind_DA_min] = min( ESS_N_sig(1,3:end));
ind_DA_min = ind_DA_min + 2;

N_ESS_DA_min = ESS_N_sig(:,ind_DA_min);

N_mean_DA_min = zeros(8,1);
N_std_DA_min = zeros(8,1);
for jj = 1:K
    xxx = eval(char(NN_all{jj}));
    xxx = xxx(ind_DA_min,:);
    N_mean_DA_min(jj) = mean(xxx);
    N_std_DA_min(jj) = std(xxx);
end


[~, ind_DA_max] = max( ESS_N_sig(1,3:end));
ind_DA_max = ind_DA_max + 2;

N_ESS_DA_max = ESS_N_sig(:,ind_DA_max);

N_mean_DA_max = zeros(8,1);
N_std_DA_max = zeros(8,1);
for jj = 1:K
    xxx = eval(char(NN_all{jj}));
    xxx = xxx(ind_DA_max,:);
    N_mean_DA_max(jj) = mean(xxx);
    N_std_DA_max(jj) = std(xxx);
end
 

% ind_sel = [5,15,25,35];
% ind_sel = [6,15,24,33];
ind_sel = 4:4:36;

N_ESS_sel =  ESS_N_sig(:,ind_sel);

N_mean_sel = zeros(8,length(ind_sel));
N_std_sel = zeros(8,length(ind_sel));
    
for jj = 1:K
    xxx = eval(char(NN_all{jj}));    
    xxx = xxx(ind_sel,:);
    N_mean_sel(jj,:) = mean(xxx,2);
    N_std_sel(jj,:) = std(xxx,0,2);    
end

% Theta
TH_ESS = ESS_TH_sig;
TH_mean = zeros(8,9);
TH_std = zeros(8,9);

for jj = 1:K
    xxx = eval(char(TH_all{jj}));
    TH_mean(jj,:) = mean(xxx);
    TH_std(jj,:) = std(xxx);
end 
 
% save
% BurnIn = 10000;
results_path = ['Results/BurnIn_',num2str(BurnIn),ext,'/'];
name = [results_path,'BKM_ALL_results_selected',ext,'.mat'];
save(name,'-regexp', '^TH', '^N_','^ind_','Time')

% table
% clear all
% ext = '_Selected';
% BurnIn = 10000;
results_path = ['Results/BurnIn_',num2str(BurnIn),ext,'/'];
name = [results_path,'BKM_ALL_results_selected',ext,'.mat'];

load(name, '-regexp','^TH', 'Time')
print_table_theta(TH_mean, TH_std, TH_ESS, Time, ext, results_path)

load(name, '-regexp','^N_', 'Time')
print_table_Na(N_mean_sel, N_std_sel, N_ESS_sel,...
    N_mean_DA_min, N_std_DA_min, N_ESS_DA_min,...
    N_mean_DA_max, N_std_DA_max, N_ESS_DA_max,...
    Time,ext,results_path)
% DATA_NAMES = {'_GSPC','_IBM','_AAPL','_MSFT','_JPM','_GE'};
DATA_NAMES = {'_GSPC','_IBM','_MSFT'};

figures_path = 'figures/Empirical/';
save_on = true;
    
for ii = [1,2] %[1,2,4]
%     y = load('../Data/Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv');
    y = load('other/Perc_Rets_GSPC_IBM_MSFT.csv');

    y = y(:,ii); % arg0 = 4 = MSFT
%     time = [1998,(2017 + 8/12)];       
    time = [2000,2018];       
%     TT = length(y);
%     T = 2000;
%     T = 2000;    
%     y = y((TT-T+1):TT)';
%     time_up = time(1) + (time(2) - time(1))*T/TT;
%     time(1) = time_up;
    data_name = DATA_NAMES{ii}; %'_MSFT';  

    Plot_data(y,time,save_on,figures_path,data_name)
end

T = 2000;
figures_path = 'figures/Simulation/';
save_on = true;
time = [1,T];
Plot_data(y,time,save_on,figures_path,'_simulation')

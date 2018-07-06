start_date = '01012000';
end_date = '31122017';
time = [2000, 2018];
tickers = {'^GSPC' 'IBM' 'MSFT'};
% tickers = {'^GSPC' 'IBM' 'AAPL' 'MSFT' 'JPM' 'GE'};
% tickers2 = {'KO', 'T', 'WMT', 'XOM'}; % NAIS: GE JP Morgan Coca-Cola AT&T Wal-Mart Exxon
% tickers3 ={'FORD','BAC','SBUX','PVH'};
% 
% GSPC = hist_stock_data(start_date, end_date, '^GSPC');
% date_gspc = GSPC.Date;
% GSPC2 = GSPC.AdjClose;
% GSPC = GSPC.Close;
% 
% AAPL = hist_stock_data(start_date, end_date, 'AAPL');
% date_aapl = AAPL.Date;
% AAPL2 = AAPL.AdjClose;
% AAPL = AAPL.Close;
% 
% IBM = hist_stock_data(start_date, end_date, 'IBM');
% date_ibm = IBM.Date;
% IBM2 = IBM.AdjClose;
% IBM = IBM.Close;

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
T = size(data_Returns,1);
csvwrite('Perc_Rets_GSPC_IBM_MSFT.csv',data_Returns);

figure(1)
subplot(3,2,1)
    plot(data_Returns(:,1))
    xlim([1,T])
subplot(3,2,2)
    histogram(data_Returns(:,1))
subplot(3,2,3)
    plot(data_Returns(:,2))
    xlim([1,T])
subplot(3,2,4)
    histogram(data_Returns(:,2))
subplot(3,2,5)
    plot(data_Returns(:,3))
    xlim([1,T])
subplot(3,2,6)
    histogram(data_Returns(:,3))

    
    
times = data.Date(2:end);
times(1) %     '1998-01-05'
times(end)  % '2017-08-30'

TT = length(times);
T = 2000;
times(TT-T+1) %  '2009-09-22'
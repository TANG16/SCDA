function fname = print_table_results_theta_and_h(data_name, method,...
    TH_mean, TH_std, TH_ESS,  H_mean, H_std, H_ESS, Time)
% data_name = '_IBM';
% T = 4527;
K = length(method);
 
ind_h_sel = 2:50:T;
ind_h_sel = ind_h_sel(10*(1:9));
h_sel = cell(1,length(ind_h_sel));
h_sel_str = '';
for ii = 1:length(ind_h_sel);
    h_sel{1,ii} = sprintf('$h_{%i}$',ind_h_sel(ii));
    h_sel_str = [h_sel_str, ' & ', sprintf('$h_{%i}$',ind_h_sel(ii))];
end
params = {'$\\mu$','$\\phi$','$\\sigma^2$', '$\\beta$', '$\\rho$'};

path = ['Results/Empirical/'];    


%% THETA
    fname = [path,'SV_results_theta',data_name,'.tex'];

    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, ['\\begin{tabular}{cc ccccc } \n']);
    fprintf(FID, [' Method & & '...
        ' $\\mu$ & ','$\\phi$ &','$\\sigma^2$','$\\beta$','\\rho',...        
        ' \\\\ \\hline  \\hline\n']);
    
    for ii = 1:K
        fprintf(FID, '\\rowcolor{LightCyan} \n');
        fprintf(FID, '%s & Mean \n', method{ii});
        for jj = 1:3
            fprintf(FID,' & %6.4f ', TH_mean(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

        fprintf(FID, ' & (Std) \n');
        for jj = 1:3
            fprintf(FID,' & (%6.4f) ', TH_std(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n'); 
        
        fprintf(FID, ' [%6.2f s] & ESS \n', Time(ii));
        for jj = 1:3
            fprintf(FID,' & %6.4f ', TH_ESS(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

%         fprintf(FID, '[%6.2f s]  & ESS/sec. \n', Time(ii));
%         for jj = 1:9
%         fprintf(FID,' & %6.4f ',TH_ESS(ii,jj)/Time(ii));
%         end
%         fprintf(FID,' \\\\  [1.3ex] \n');        
    end
    fprintf(FID,'  \\hline \n');


    fprintf(FID, ['\\multicolumn{7}{p{6cm}}{\\footnotesize{ESS: at lag equal to the ' ...
    'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{7}{p{6cm}}{\\footnotesize{Computing times '...
        ' (in seconds) in square brackets.}}  \\\\ \n']);

    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{SV model: posterior means, standard deviations and '...
        'effective sample sizes (ESS) of the model parameters ' ...
        'for $M=10000$ posterior draws after a burn-in of $10000$.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:SVML_results_theta',data_name,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fprintf(FID, '} \\normalsize');
    fclose(FID);
    
    
    
%% H
    fname = [path,'SV_results_H',data_name,'.tex'];

    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, ['\\begin{tabular}{cc ', repelem('c',length(ind_h_sel)),'} \n']);
    fprintf(FID, [' Method & ', h_sel_str, ... 
        ' \\\\ \\hline  \\hline\n']);
    
    for ii = 1:K
        fprintf(FID, '\\rowcolor{LightCyan} \n');
        fprintf(FID, '%s & Mean \n', method{ii});
        for jj = 1:length(ind_h_sel)
            fprintf(FID,' & %6.4f ', H_mean(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

        fprintf(FID, ' & (Std) \n');
        for jj = 1:length(ind_h_sel)
            fprintf(FID,' & (%6.4f) ', H_std(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n'); 
        
        fprintf(FID, ' [%6.2f s] & ESS \n', Time(ii));
        for jj = 1:length(ind_h_sel)
            fprintf(FID,' & %6.4f ', H_ESS(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

%         fprintf(FID, '[%6.2f s]  & ESS/sec. \n', Time(ii));
%         for jj = 1:9
%         fprintf(FID,' & %6.4f ',TH_ESS(ii,jj)/Time(ii));
%         end
%         fprintf(FID,' \\\\  [1.3ex] \n');        
    end
    fprintf(FID,'   \\hline \n');


    fprintf(FID, ['\\multicolumn{5}{p{8cm}}{\\footnotesize{ESS: at lag equal to the ' ...
    'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{5}{p{8cm}}{\\footnotesize{Computing times '...
        ' (in seconds) in square brackets.}}  \\\\ \n']);

    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{SV model: posterior means, standard deviations and '...
        'effective sample sizes (ESS) of the latent volatilities ' ...
        'for $M=10000$ posterior draws after a burn-in of $10000$.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:SV_results_H',data_name,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fprintf(FID, '} \\normalsize');
    fclose(FID);
    
end
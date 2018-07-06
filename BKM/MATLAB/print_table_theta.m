function fname = print_table_theta(TH_mean, TH_std, TH_ESS, Time,ext,path)

    fname = [path,'BKM_theta',ext,'.tex'];
    
%     params = {'$\\alpha_{1}$', '$\\alpha_{a}$', '$\\alpha_{\\rho}$', '$\\alpha_{\\lambda}$', ...
%     '$\\beta_{1}$', '$\\beta_{a}$', '$\\beta_{\\rho}$', '$\\beta_{\\lambda}$',...
%     '$\\sigma^{2}_{y}$'};

    method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30', 'Exact'};
    K = length(method);
    
    FID = fopen(fname, 'w+');
    fprintf(FID, '{\\setlength{\\tabcolsep}{2pt}\n');

    fprintf(FID, '{\\footnotesize \n');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
%     fprintf(FID, '\\hspace*{-2.5cm} \n');
    fprintf(FID, '\\begin{tabular}{cc ccc ccc ccc} \n');
    fprintf(FID, '\\hline \n');
    fprintf(FID, [' Method & & '...
        ' $\\alpha_{1}$ & ', ' $\\alpha_{a}$ & ', ' $\\alpha_{\\rho}$ & ', ' $\\alpha_{\\lambda}$ & ', ...
        ' $\\beta_{1}$ & ', ' $\\beta_{a}$ & ', ' $\\beta_{\\rho}$ & ', ' $\\beta_{\\lambda}$ & ',...
        ' $\\sigma^{2}_{y}$ '...
        ' \\\\ \\hline  \\hline\n']);
    
    for ii = 1:K
        fprintf(FID, '\\rowcolor{LightCyan} \n');
        fprintf(FID, '%s & Mean \n', method{ii});
        for jj = 1:9
        fprintf(FID,' & %6.4f ', TH_mean(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

        fprintf(FID, ' & (Std) \n');
        for jj = 1:9
        fprintf(FID,' & (%6.4f) ', TH_std(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n'); 
        
        fprintf(FID, ' & ESS \n');
        for jj = 1:9
        fprintf(FID,' & %6.4f ', TH_ESS(ii,jj));
        end
        fprintf(FID,' \\\\  [0.75ex] \n');

        fprintf(FID, '[%6.2f s]  & ESS/sec. \n', Time(ii));
        for jj = 1:9
        fprintf(FID,' & %6.4f ',TH_ESS(ii,jj)/Time(ii));
        end
        fprintf(FID,' \\\\  [1.3ex] \n');        
    end
    fprintf(FID,' \\\\  \\hline \n');

    fprintf(FID, ['\\multicolumn{11}{p{11cm}}{\\footnotesize{ESS: at lag equal to the ' ...
    'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{11}{p{11cm}}{\\footnotesize{Computing times '...
        ' (in seconds) in square brackets.}}  \\\\ \n']);

    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Posterior means, standard deviations and '...
        'effective sample sizes (ESS) of the model parameters ' ...
        'for $M=10000$ posterior draws after a burn-in of $10000$ for the lapwings data.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:BKM_theta',ext,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fprintf(FID, '} \\normalsize');
    fclose(FID);
end
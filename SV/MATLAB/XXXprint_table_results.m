function fname = print_table_results(TH_mean, TH_std, TH_ESS, TH_AR,...
    HH_mean, HH_std, HH_ESS, HH_AR,...
    Time,data,data_name,path)


    params = {'$\mu$','$\phi$','$\sigma^2$'};
    ind_sel = {'$h_{200}$','$h_{600}$','$h_{1000}$','$h_{1400}$','$h_{1800}$'};
    
    fname = [path,'/SV_results',data,'.tex'];
%     lags = {'sig','40', '100','1000'};
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{ccc cc cc} \n');
    fprintf(FID, '\\hline \n');
%     fprintf(FID, [' & ESS lag&& DA KSC & DA RW & HMM fix & HMM adapt \\\\ \\hline  \\hline\n']);
    fprintf(FID, [' & && DA KSC & DA RW & HMM fix & HMM adapt \\\\ \\hline  \\hline\n']);
    
    fprintf(FID, ' \\multicolumn{2}{c}{time (s)}&');
    for ii = 1:4
        fprintf(FID,' & %6.4f ', Time(ii)); 
    end
    fprintf(FID,' \\\\  \\hline \n');

    fprintf(FID, ['\\multicolumn{7}{c}{$\\theta$} \\\\ \\hline \n']);  
    for pp = 1:3
        fprintf(FID,'\\multirow{4}{*}{%s}  ',params{pp});

        fprintf(FID,' & mean  &  ');        
        for ii = 1:4
            fprintf(FID,' & %6.4f ', TH_mean(ii,pp));            
        end
        fprintf(FID,' \\\\ [0.75ex]\n');  
        
        fprintf(FID,' & std  &  ');        
        for ii = 1:4
            fprintf(FID,' & (%6.4f) ', TH_std(ii,pp));            
        end
        fprintf(FID,' \\\\ [0.75ex]\n');  
        
        fprintf(FID,' & ESS  &  ');        
        for ii = 1:4
            fprintf(FID,' & %6.4f ', TH_ESS(ii,pp));            
        end
        fprintf(FID,' \\\\ [0.75ex]\n');        
        
        fprintf(FID,' & AR && -- ');
        for ii = 2:4
            fprintf(FID,' & %6.4f ', TH_AR(ii,pp)); 
        end
        fprintf(FID,' \\\\ [1.3ex] \n');        
    end
    fprintf(FID, '\\hline \n');

    fprintf(FID, ['\\multicolumn{7}{c}{$ \\bm{h} $} \\\\ \\hline \n']);  
    for pp = 1:5
        fprintf(FID,'\\multirow{4}{*}{%s}  ',ind_sel{pp});

        fprintf(FID,' & mean &  ');        
        for ii = 1:4
            fprintf(FID,' & %6.4f ', HH_mean(ii,pp));            
        end
        fprintf(FID,' \\\\ [0.75ex]\n');  
        
        fprintf(FID,' & std &  ');        
        for ii = 1:4
            fprintf(FID,' & (%6.4f) ', HH_std(ii,pp));            
        end
        fprintf(FID,' \\\\ [0.75ex]\n');  
        
        fprintf(FID,' & ESS  &  ');        
        for ii = 1:4
            fprintf(FID,' & %6.4f ', HH_ESS(ii,pp));            
        end       
        fprintf(FID,' \\\\ [1.3ex] \n');        
    end
    
    fprintf(FID,' & AR & & -- ');
    for ii = 2:4
        fprintf(FID,' & %6.4f ', HH_AR(ii)); 
    end
    fprintf(FID, '\\\\ \\hline \n');          
        
    fprintf(FID, '\\hline \n');
 
    fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{ESS: effective sample size '...
        'at lag equal to the ' ...
        'lowest order at which sample autocorrelation is not significant.}}  \\\\ \n']);
    fprintf(FID, ['\\multicolumn{7}{p{11cm}}{\\footnotesize{AR: acceptance rate ' ...
    'of the MH RW algorithm. (for $h_{t}$ average over imputations)}}  \\\\ \n']);
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Posterior means, standard deviations, ',...
        'effective sample sizes (ESS) ',...
        ' for $M=10000$ posterior draws after a burn-in of $10000$ ' ,...
        'for $T=2000$ observations of \\textbf{',data_name ,'} data.}\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:SV_results',data,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end
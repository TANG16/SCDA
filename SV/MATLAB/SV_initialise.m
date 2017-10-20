function [theta, delta] = SV_initialise(arg0, method)

    theta_init = [0, 0.95, 0.1^2];
    theta_init2 = [0, 0.95, 0.16^2];
    theta_init3 = [0, 0.95, 0.4^2];
    delta.h = 0.1;
    delta.t = [0.4, 0.01, 0.004];

    switch arg0
        case 0 % simulation
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
    %         if arg1       
                    delta.h = 0.2;               
                    theta = theta_init;   
    %         end
                case 2  % FULL DA - RANDOM WALK UPDATES  
    %         if (arg2 == 1)  
    %             delta.h = 0.15;
    %             theta = theta_init2;        
    %         elseif (arg2 == 2)
                    delta.h = 0.15;
                    theta = theta_init2;
    %         end
                case 3 % SEMI DA:
    %         if arg3
                    delta.t = [0.3, 0.01, 0.004];
    %             if (N_bin >= 30)
                    delta.h = 0.5;
    %             else
    %                 delta.h = 0.1;            
    %             end
                    theta = theta_init2;
    %         end
                case 4 % SEMI DA shift: 
    %         if arg4
                    delta.t = [0.4, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init3;
    %         end
                case 5 % SEMI DA efficient implementation:
    %         if ((nargin == 6) || ((nargin == 7) && (arg5 == 1)))
    %             if (N_bin >= 30)
                    delta.t = [0.3, 0.01, 0.004]; 
                    delta.h = 0.5;
    %             else
    %                 delta.t = [0.03, 0.01, 0.004]; 
    %                 delta.h = 0.1;            
    %             end
                    theta = theta_init2;     
    %         end
                otherwise % adaptive
                    delta.t = [0.1, 0.01, 0.001];
                    delta.h = 0.1;
                    theta = theta_init;
            end
        case 1 % GSPC
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
                    delta.h = 0.55; % A_H 0.38 % 0.45 - A_H 0.43 %%% 0.45 - A_H 0.34; %0.5 - A_H 0.30; % 0.2 - A_H_DA 0.62             
                    theta = theta_init2; %theta_init   
                case 2  % FULL DA - RANDOM WALK UPDATES  
                    delta.h = 0.65; % A_H 0.36; % 0.55 - A_H 0.41 %%% 0.55 - A_H - 0.33; 0.5 - A_H 0.39  % 0.25 - A_H 0.62; % 0.15 - A_H 0.75
% %                     delta.t = [0.4, 0.01, 0.004]; %     0.3332    0.4781    0.4454
% %                     delta.t = [0.4, 0.015, 0.005];  %     0.3249    0.3628    0.3658
%                     delta.t = [0.4, 0.015, 0.005]; %     0.3161    0.3802    0.4271
                    delta.t = [0.4, 0.015, 0.007]; %    0.3153    0.3839    0.3341
                    theta = theta_init2;
                case 3 % SEMI DA:
                    delta.t = [0.3, 0.01, 0.004];
                    delta.h = 0.5;
                    theta = theta_init2;
                case 4 % SEMI DA shift: 
                    delta.t = [0.4, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init3;
                case 5 % SEMI DA efficient implementation:
%                     delta.t = [0.3, 0.01, 0.004];  %     0.3752    0.5037    0.6074
%                     delta.t = [0.35, 0.02, 0.007]; %     0.3345    0.3172    0.4482
                    delta.t = [0.35, 0.02, 0.010]; %     0.3280    0.3119    0.3492
                    delta.h = 0.86; % 0.36 %0.83 - A_H 0.38; 0.83; % 0.8 - A_H 0.39; %0.79 - A_H 0.39; % 0.77 - A_H 0.41; %0.8 - A_H 0.05 % 0.75 - A_H 0.42; % 0.6 - A_H 0.49; % 0.5 - A_H 0.55
                    theta = theta_init2;     
                otherwise % adaptive
                    delta.t = [0.1, 0.01, 0.001];
                    delta.h = 0.2;
                    theta = theta_init2; % [0, 0.95, 0.16^2];
                    theta(3) =  1;
            end
            
        case 4 % MSFT
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
                    delta.h = 0.2;               
                    theta = theta_init;   
                case 2  % FULL DA - RANDOM WALK UPDATES  
                    delta.h = 0.15;
                    theta = theta_init2;
                case 3 % SEMI DA:
                    delta.t = [0.3, 0.01, 0.004];
                    delta.h = 0.5;
                    theta = theta_init2;
                case 4 % SEMI DA shift: 
                    delta.t = [0.4, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init3;
                case 5 % SEMI DA efficient implementation:
                    delta.t = [0.3, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init2;     
                otherwise % adaptive
                    delta.t = [0.1, 0.01, 0.001];
                    delta.h = 0.1;
%                     theta = [0.1, 0.98, 0.1];
                    theta = [ 0.3650    0.8717    0.1586];
            end
            
        otherwise
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
                    delta.h = 0.2;               
                    theta = theta_init;   
                case 2  % FULL DA - RANDOM WALK UPDATES  
                    delta.h = 0.15;
                    theta = theta_init2;
                case 3 % SEMI DA:
                    delta.t = [0.3, 0.01, 0.004];
                    delta.h = 0.5;
                    theta = theta_init2;
                case 4 % SEMI DA shift: 
                    delta.t = [0.4, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init3;
                case 5 % SEMI DA efficient implementation:
                    delta.t = [0.3, 0.01, 0.004]; 
                    delta.h = 0.5;
                    theta = theta_init2;     
                otherwise % adaptive
                    delta.t = [0.1, 0.01, 0.001];
                    delta.h = 0.1;
                    theta = theta_init;
            end 
    end
end
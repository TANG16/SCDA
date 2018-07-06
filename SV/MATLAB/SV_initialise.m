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
%                     delta.h = 0.2;     % 0.5825          
%                     delta.h = 0.4;     % 0.3896          
                    delta.h = 0.5;     % 0.3194        
                    theta = theta_init2;   %  
    %         end
                case 2  % FULL DA - RANDOM WALK UPDATES  
%                     delta.h = 0.15; %     0.6625
                    delta.h = 0.4; %     0.0.3696
%                     delta.t = [0.4, 0.01, 0.004]; %    0.4667    0.3524    0.2516
                    delta.t = [0.6, 0.01, 0.003]; %      0.3603    0.3472    0.3157
                    theta = theta_init2;
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
%                     delta.t = [0.3, 0.01, 0.004]; %    0.3764    0.3493    0.3383
%                     delta.t = [0.35, 0.01, 0.004]; %    0.0036    0.1802    0.0149
%                     delta.t = [0.25, 0.01, 0.004]; %      0.4529    0.3641    0.3481
%                     delta.t = [0.28, 0.01, 0.004]; %        0.4020    0.3550    0.3406
%                     delta.t = [0.3, 0.01, 0.004]; %       0.0041    0.1043    0.0180                     
%                     delta.t = [0.29, 0.01, 0.004]; %      0.0040    0.1652    0.0146
%                     delta.h = 0.5; % 0.4038
%                     delta.h = 0.6; % 0.0633
%                     delta.h = 0.55; %    0.0778
%                     delta.h = 0.5; %      0.4107
%                     delta.h = 0.55; %      0.3674 %    0.0683
%                     delta.h = 0.56; %         0.3702
%                     delta.h = 0.57; %   0.0761   
                    delta.t = [0.3, 0.01, 0.004]; %        0.0041    0.1043    0.0180
                    delta.h = 0.5; 
                    theta = theta_init2;     
                 otherwise % adaptive
%                     delta.t = [0.1, 0.01, 0.001]; %    0.7386    0.3556    0.7117
%                     delta.t = [0.4, 0.01, 0.004]; %    0.3800    0.3524    0.3365
%                     delta.t = [0.45, 0.01, 0.010];  %    0.3448    0.3538    0.1585
                    delta.t = [0.45, 0.01, 0.004];  %   0.3509    0.3503    0.3396
%                     delta.h = 0.1;  %  0.8340
%                     delta.h = 0.51;  %   0.3977      
%                     delta.h = 0.55;  %    0.3792
                    delta.h = 0.6;  %  0.3510
                                  
                    theta = theta_init2;
            end
        case 1 % GSPC
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
%                     delta.h = 0.55; %0.4063
%                     delta.h = 0.6; %0.3793           
                    delta.h = 0.7; %0.3336           
                    theta = theta_init2; %theta_init   
                case 2  % FULL DA - RANDOM WALK UPDATES  
%                     delta.h = 0.65;  %    0.3759
                    delta.h = 0.7;  %      0.3478
%                     delta.t = [0.4, 0.015, 0.007]; %      0.3104    0.3905    0.3601
                    delta.t = [0.4, 0.018, 0.007]; %        0.3106    0.3381    0.3483
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
%                     delta.t = [0.25, 0.02, 0.008];  %         0.4252    0.3196    0.4291
%                     delta.h = 0.7; %      0.4511
%                     delta.t = [0.25, 0.02, 0.008];  %            0.4275    0.3147    0.4102
%                     delta.h = 0.8; %      0.3995
                    delta.t = [0.25, 0.02, 0.008];  %            0.4275    0.3147    0.4102
                    delta.h = 0.8; %      0.3995                    
                    theta = theta_init2;     
                otherwise % adaptive
%                     delta.t = [0.1, 0.01, 0.001];     0.5614    0.7170    0.9383
                    delta.t = [0.35, 0.02, 0.010];   %         0.3360    0.3114    0.3488             
%                     delta.h = 0.2; %     0.9124
%                     delta.h = 0.9; %     0.3629
                    delta.h = 1.0; %   0.3392
                    theta = theta_init2; % [0, 0.95, 0.16^2];
%                     theta(3) =  1;
            end
            
        case 4 % MSFT
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)                  
%                     delta.h = 0.2; % 0.8326               
%                     delta.h = 0.8; % 0.4633
%                     delta.h = 1.0; % 0.4068                   
%                     delta.h = 1.2; % 0.3608                   
                    delta.h = 1.3; % 0.3162                   
                    theta = theta_init2;   
                case 2  % FULL DA - RANDOM WALK UPDATES 
%                     delta.t = [0.6, 0.01, 0.003]; %    0.0958    0.6726    0.8120
%                     delta.h = 0.15; % 0.8749
%                     delta.t = [0.26, 0.04, 0.01]; %      0.2235    0.2991    0.5354                   
%                     delta.h = 0.6;    % 0. 5691
%                     delta.t = [0.2, 0.03, 0.03];     %     0.2868    0.3728    0.2393               
%                     delta.h = 0.9;  % 0.4339                      
%                     delta.t = [0.16, 0.03, 0.025];     %        0.3407    0.3685    0.2690             
%                     delta.h = 1.2;  % 0.3419 
%                     delta.t = [0.16, 0.03, 0.02];     %        0.3751    0.2921    0.1753
%                     delta.h = 1.2;  % 0.2253
%                     delta.t = [0.16, 0.03, 0.02];     %     0.3420    0.3725    0.3254    
%                     delta.h = 1.1;  % 0.3419
                    delta.t = [0.16, 0.03, 0.02];     %      0.3469    0.3558    0.3023   
                    delta.h = 1.1; % 0.3448
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
%                     delta.t = [0.10, 0.04, 0.020]; %     0     
%                     delta.h = 1.0; %  0
%                     delta.t = [0.10, 0.04, 0.020]; %   0.4652    0.3155    0.4504    
%                     delta.h = 0.7; %  0.6063                    
%                     delta.t = [0.12, 0.04, 0.025]; %   0.4113    0.3189    0.3913   
%                     delta.h = 0.8;   % 0.5637           
%                     delta.t = [0.13, 0.04, 0.025]; %  0  
%                     delta.h = 0.9;   % 0                   
%                     delta.t = [0.13, 0.04, 0.025]; %      0.3973    0.3155    0.3750
%                     delta.h = 0.85;  % 0.53          
%                     delta.t = [0.10, 0.03, 0.020]; %     0.4606    0.4054    0.4735 
%                     delta.h = 1.0;
%                     delta.t = [0.10, 0.03, 0.020]; %     0.4606    0.4054    0.4735 
%                     delta.h = 1.1; 
%                     delta.t = [0.050, 0.02, 0.015]; %        0.6443    0.5231    0.5452
%                     delta.h = 0.7;           % 0.6135         
                    delta.t = [0.050, 0.02, 0.015]; %        
                    delta.h = 0.9;                    
                    theta = theta_init2;     
                otherwise % adaptive
%                     delta.t = [0.3, 0.01, 0.004];  %    0.2028    0.6820    0.8090
%                     delta.h = 0.7; %    0.5986
%                     delta.t = [0.25, 0.03, 0.01];  %     0.2353    0.3910    0.6311   
%                     delta.h = 1.1;    %     0.4631
%                     delta.t = [0.15, 0.04, 0.03];  %     0.3634    0.3157    0.3394 
%                     delta.h = 1.5;   %     0.3680              
                    delta.t = [0.17, 0.04, 0.03];  %      0.3351    0.3085    0.3116
                    delta.h = 1.7;  %0.3232
                    %                     theta = [0.1, 0.98, 0.1];
                    theta = theta_init2;
            end
            
        otherwise
            switch method 
                case 1  % FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
%                     delta.h = 0.2;    % 0.8180
%                     delta.h = 0.6;    % 0.4695  
%                     delta.h = 1.0;      %    0.1287     
%                     delta.h = 0.8;   % 0.3535  --> 0.4572
%                     delta.h = 1.0; % 0.3857
%                     delta.h = 1.1; % 0.3622
                    delta.h = 1.2; % 0.3404
                    theta = theta_init2;   
                case 2  % FULL DA - RANDOM WALK UPDATES  
%                     delta.t = [0.4, 0.01, 0.004];  %  0.1557    0.6621    0.7407               
%                     delta.t = [0.25, 0.03, 0.010]; %  0.2335    0.3439    0.4466
%                     delta.t = [0.2, 0.03, 0.015]; %  0.2891    0.3711    0.4190
%                     delta.t = [0.15, 0.03, 0.02]; %  0.3660    0.3331    0.2610
%                     delta.t = [0.16, 0.03, 0.015]; %      0.3520    0.3185    0.2921 
%                                                    %   -->  0.1689    0.2647    0.0258
%                     delta.t = [0.16, 0.03, 0.015]; %     0.3397    0.3262    0.2883
%                     delta.t = [0.16, 0.03, 0.014]; %       0.3569    0.3624    0.4223
                    delta.t = [0.16, 0.03, 0.015]; %   0.3397    0.3262    0.2883
%                     delta.h = 0.15; %     0.8686
%                     delta.h = 0.6; % 0.5471
%                     delta.h = 1.0; %     0.3557          0.4066           
%                     delta.h = 1.2; %    0.2949        0.2581   
%                     delta.h = 1.0; %    0.3615   0.3138 --> 0.0831
                    delta.h = 0.9; % 0.3372  
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
% %                     delta.t = [0.2, 0.04, 0.015];  %    0.2926    0.3107    0.5189
% %                     delta.t = [0.18, 0.04, 0.015];  %       0.3241    0.3100    0.5153
% %                     delta.t = [0.18, 0.04, 0.02];  %       0.0058    0.1209    0.0025  
% %                     delta.t = [0.18, 0.04, 0.018];  %      0.0060    0.1086    0.0031 
% %                     delta.t = [0.18, 0.04, 0.017]; %    0.3259    0.3088    0.4770
% %                     delta.t = [0.18, 0.04, 0.017]; %   0.0051    0.1189    0.0022
% %                     delta.t = [0.18, 0.04, 0.03]; %        0.3240    0.3044    0.3234
% %                     delta.h = 0.9; %    0.5251 
% %                     delta.h = 1.1; %    0.4560 
% %                     delta.h = 1.5; %       0.0248  
% %                     delta.h = 1.3;   % 0.0279                 
% %                     delta.h = 1.1;   %               0.4550       
% %                     delta.h = 1.2;   %               0.4304
% %                     delta.h = 1.55;   %            0.3524
% %                     delta.h = 1.65;   %            0.3350
% %                     delta.h = 1.75;   %            0.31750
%                     delta.t = [0.18, 0.04, 0.03]; %     0.0040    0.1078    0.0007
%                     delta.h = 1.75; % 0.0153
                    delta.t = [0.17, 0.04, 0.03]; %        0.3413    0.3069    0.3220
                    delta.h = 1.5; % 0.3584
                    theta = theta_init2;     
                otherwise % adaptive             
%                     delta.t = [0.17, 0.04, 0.03];  % 0.5188    0.6304    0.9225
                    delta.t = [0.17, 0.04, 0.03];   %     0.3420    0.2995    0.3219
                    delta.h = 1.7;   % 0.3274
                    theta = theta_init2;
            end
    end
end
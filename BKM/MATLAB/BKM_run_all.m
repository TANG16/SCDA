addpath(genpath('../MATLAB'));

long = true;

if ~long 
    M = 10000; %50000; %10000;
    save_on = false; %true;

    % for iter = 4:4    
        for BurnIn = 10000; %20000 % [0, 1000, 5000, 10000]
    %         if (iter > 4    )
                if save_on
                    oldname = ['Results/BurnIn_',num2str(BurnIn)];
                    mkdir(oldname)
                end
                Results_DA = BKM_try_DA(M, BurnIn, save_on);
                Results_Exact = BKM_try_HMM(M, BurnIn, save_on);
    %         end

            ii = 0;
            for B = [10,20,30] % [30]; %[10,20,30]
    %             cd Adaptive
                ii = ii + 1;

                OK = false;
                while (OK == false)
                    try
                        Results_Adapt{ii} = BKM_try_HMM_adapt(B, M, BurnIn, save_on);
                        OK = true;
                    catch
        %                 Results_Adapt{ii} = BKM_try_HMM_adapt(B, M, BurnIn, save_on);
                    end
                end
    %             cd ../Bins
                OK = false;
                while (OK == false)
                    try 
                        Results_Bin{ii} = BKM_try_HMM_bin(B, M, BurnIn, save_on);
                        OK = true;
                    catch
        %                 Results_Bin{ii} = BKM_try_HMM_bin(B, M, BurnIn, save_on);
                    end
                end
    %             cd ..
            end

            if save_on
                newname = ['Results/BurnIn_',num2str(BurnIn),'_',num2str(iter)];
                movefile(oldname,newname)
            end
        end
    % end
else
    M = 100000; %50000; %10000;
    save_on = true;

    BurnIn = 10000;
    BKM_try_DA(M, BurnIn, save_on);
%     BKM_try_HMM(M, BurnIn, save_on);
%     
%     BKM_try_HMM_bin(30, M, BurnIn, save_on);
%     BKM_try_HMM_bin(20, M, BurnIn, save_on);
%     BKM_try_HMM_bin(10, M, BurnIn, save_on);
    
    
%     BKM_try_HMM_adapt(30, M, BurnIn, save_on);
%     BKM_try_HMM_adapt(20, M, BurnIn, save_on);
%     BKM_try_HMM_adapt(10, M, BurnIn, save_on);
    
end
% BKM_try_DA(M,5000)

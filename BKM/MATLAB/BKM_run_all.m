M = 10000;

for BurnIn = 10000 % [0, 1000, 5000, 10000]
    BKM_try_DA(M,BurnIn)
    BKM_try_HMM(M,BurnIn)
    for B = [10,20,30]
        cd Adaptive
        try
            BKM_try_HMM_adapt(B, M, BurnIn)
        catch
            BKM_try_HMM_adapt(B, M, BurnIn)
        end
        cd ../Bins
        try 
            BKM_try_HMM_bin(B, M, BurnIn)
        catch
            BKM_try_HMM_bin(B, M, BurnIn)
        end
        cd ..
    end
end

% BKM_try_DA(M,5000)

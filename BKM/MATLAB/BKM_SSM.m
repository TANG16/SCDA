function KFS_params = BKM_SSM(phi1,phia,rho,sigy)

    n = max(size(phi1));
    
    KFS_params.c = 0;
    KFS_params.Z = [0,1];
    KFS_params.H = sigy;

    KFS_params.d = [0; 0];
    T = zeros(2,2,n);
    T(1,2,:) = rho.*phi1;
    T(2,1,:) = phia;
    T(2,2,:) = phia;
    KFS_params.T = T;

    Q = zeros(2,2,n);
    Q(1,1,:) = rho.*phi1;
    Q(2,2,:) = (1-phia).*phia;
    KFS_params.Q = Q;

    KFS_params.R = eye(2,2);
    KFS_params.R1 = [0,1;0,0];
    KFS_params.R2 = [0,0;1,1];
    
    KFS_params.a1 = [200; 1000];
    KFS_params.P1 = 1000*eye(2,2);
end
%************************** NOTES **************************
% 1. Initialization points matter; 
%    i. RSI is small or big -- powerMatInitOne(average) & powerMatInitTwo
%    (waterfilling of relay link) & powerMatInitThree (waterfilling of direct link)
%    ii. Time-dimension allocations are the same -- powerMatInitFour(time half-duplex)
%***********************************************************
%% 
clear;

%%
N       = 4; % total number of subcarriers
Delta_P = 30;% dB. The average power per-subcarrier
delta_P = 10.^(Delta_P/10);
P_S     = 0.5*delta_P*N; % Source power
P_R     = 0.5*delta_P*N; % Destination power
omega_sr = 1;
omega_rr = 1;
omega_rd = 1;
omega_sd = sqrt(0.01);
Rho_dB   = -20:2.5:10;
Rho     = 10.^(Rho_dB/10);
T       = 100; % channel realization number

rateScheme_2 = zeros(length(Rho),1);
rate = zeros(length(Rho),1);
iterMat2 = zeros(length(Rho),T);
iterMat3 = zeros(length(Rho),T);
iterMat4 = zeros(length(Rho),T);

% rate_data_scheme_2 = zeros(length(Rho),1);
% save('rate_data.mat','rate_data_scheme_2','-append');

tic
tt = 1;
while(tt<=T)
    tt
    Hsr     = omega_sr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
    Hrr     = omega_rr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
    Hrd     = omega_rd/sqrt(2)*(randn(N,1)+1i*randn(N,1));
    Hsd     = omega_sd/sqrt(2)*(randn(N,1)+1i*randn(N,1));
    
    gammaSRT = abs(Hsr).^2;
    gammaRDT = abs(Hrd).^2;
    gammaSDT = abs(Hsd).^2;
    
    
    powerMatInit3 = powerMatInitThree(P_S,N,gammaSDT);
    powerMatInit4 = powerMatInitFour(P_S,P_R,N,gammaSRT,gammaRDT);
    timeRatioInit = 0.5;
    
    
    gammaRRT = zeros(N,length(Rho));
    feasId = zeros(length(Rho),3);
    
%     for kk = 1:length(Rho)
%         
%         powerMatInit2 = powerMatInitTwo(P_S,P_R,N,gammaSRT,gammaRDT,gammaRRT(:,kk),gammaSDT);
%         feasId(kk,1) = scheme2_SCP_CheckFeasibility(timeRatioInit,powerMatInit2,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
%         feasId(kk,2) = scheme2_SCP_CheckFeasibility(timeRatioInit,powerMatInit3,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
%         feasId(kk,3) = scheme2_SCP_CheckFeasibility(timeRatioInit,powerMatInit4,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
%     end
%     
%     if sum(sum(feasId)) ~= 3*length(Rho)
%         continue;
%     end
    
    iterRow2 = zeros(length(Rho),1);
    iterRow3 = zeros(length(Rho),1);
    iterRow4 = zeros(length(Rho),1);
    
    parfor kk = 1:length(Rho)
        
        gammaRRT(:,kk) = Rho(kk)*abs(Hrr).^2;
        powerMatInit2 = powerMatInitTwo(P_S,P_R,N,gammaSRT,gammaRDT,gammaRRT(:,kk),gammaSDT);
        
        % initial point 2 
        [rate2,timeRatio2,powerMat2,iter2] = alterOptim2(timeRatioInit, powerMatInit2,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
        
        % initial point 3
        [rate3,timeRatio3,powerMat3,iter3] = alterOptim2(timeRatioInit, powerMatInit3,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
        
        % initial point 4
        [rate4,timeRatio4,powerMat4,iter4] = alterOptim2(timeRatioInit, powerMatInit4,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT(:,kk));
        
        rate(kk) = max(max(rate2,rate3),rate4);
%         rate = rate2;
        
        iterRow2(kk) = rate2;
        iterRow3(kk) = rate3;
        iterRow4(kk) = rate4;
        
    end
    rateScheme_2 = rate + rateScheme_2; 
    
    iterMat2(:,tt) = iterRow2;
    iterMat3(:,tt) = iterRow3;
    iterMat4(:,tt) = iterRow4;
    
    tt = tt + 1;

end
rateAve_2 = rateScheme_2./T;
toc
    












































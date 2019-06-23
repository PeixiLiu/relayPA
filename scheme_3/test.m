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
Rho_dB   = -20;
Rho     = 10.^(Rho_dB/10);
T       = 50; % channel realization number

rateAve = zeros(length(Rho),1);
for tt = 22:22
    tt
    tic
%     Hsr     = omega_sr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
%     Hrr     = omega_rr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
%     Hrd     = omega_rd/sqrt(2)*(randn(N,1)+1i*randn(N,1));
%     Hsd     = omega_sd/sqrt(2)*(randn(N,1)+1i*randn(N,1));    
%     gammaSR = abs(Hsr).^2;
%     gammaRD = abs(Hrd).^2;
%     gammaSD = abs(Hsd).^2;
    for kk = 1:length(Rho)
        rate = zeros(length(Rho),1);
%         gammaRR = Rho(kk)*abs(Hrr).^2;
%         gammaSRT = zeros(2*N,1);
%         gammaRDT = zeros(2*N,1);
%         gammaSDT = zeros(2*N,1);
%         gammaRRT = zeros(2*N,1);
%         for ii = 1:2*N
%             gammaSRT(ii) = gammaSR(round(ii/2));
%             gammaRDT(ii) = gammaRD(round(ii/2));
%             gammaSDT(ii) = gammaSD(round(ii/2));
%             gammaRRT(ii) = gammaRR(round(ii/2));
%         end

        gammaSR = [1.187777377977676;0.840742630566682;1.156639843516962;0.629259961690952];
        gammaSD = [0.005826728184178;0.055808449623849;0.006056897820238;0.017477799532495];
        gammaRR = [0.001811047243540;7.173204461264228e-04;0.014475147364983;0.006857827406183];
        gammaRD = [0.139400440632621;0.193257182874610;0.462797174154918;0.989827952891790];
        
        gammaSRT = [1.187777377977676;1.187777377977676;0.840742630566682;0.840742630566682;1.156639843516962;1.156639843516962;0.629259961690952;0.629259961690952];
        gammaSDT = [0.005826728184178;0.005826728184178;0.055808449623849;0.055808449623849;0.006056897820238;0.006056897820238;0.017477799532495;0.017477799532495];
        gammaRRT = [0.001811047243540;0.001811047243540;7.173204461264228e-04;7.173204461264228e-04;0.014475147364983;0.014475147364983;0.006857827406183;0.006857827406183];
        gammaRDT = [0.139400440632621;0.139400440632621;0.193257182874610;0.193257182874610;0.462797174154918;0.462797174154918;0.989827952891790;0.989827952891790];

        powerMatInit = powerMatInitTwo(P_S,P_R,N,gammaSR,gammaRD);
        %         powerMatInit3 = powerMatInitThree(P_S,N,gammaSD);
%         powerMatInit4 = powerMatInitFour(P_S,P_R,N,gammaSR,gammaRD);
%         [powerMat2,cvx_optval2] = scheme1_SCP(powerMatInit2,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
        cvx_optvalInit = 0;
        optvalEor = 10;
        while(abs(optvalEor) > 1e-3)
               cvx_begin
                   variable powerMat(3,2*N)
                   expression obj(2*N+1)
                   expression cons(2*N+1)
                   obj(1) = 0;
                   cons(1) = 0;
                   for ii = 1:2*N
                      obj(ii+1) =  -rel_entr(1,1 + powerMatInit(3,ii)*gammaRRT(ii))...
                          - rel_entr(1,1 + powerMatInit(1,ii)*gammaSDT(ii))...
                          + rel_entr(1,1 + powerMat(3,ii)*gammaRRT(ii) + powerMat(1,ii)*gammaSRT(ii))...
                          + rel_entr(1,1 + powerMat(1,ii)*gammaSDT(ii) + powerMat(2,ii)*gammaSDT(ii))...
                          + gammaRRT(ii)/(1 + powerMatInit(3,ii)*gammaRRT(ii))*(powerMat(3,ii) - powerMatInit(3,ii))...
                          + gammaSDT(ii)/(1 + powerMatInit(1,ii)*gammaSDT(ii))*(powerMat(1,ii) - powerMatInit(1,ii))...
                          + obj(ii);
                      cons(ii+1) =  rel_entr(1,1 + powerMat(3,ii)*gammaRRT(ii))...
                          + rel_entr(1,1 + powerMat(3,ii)*gammaRDT(ii) + powerMat(1,ii)*gammaSDT(ii) + powerMat(2,ii)*gammaSDT(ii))...
                          - rel_entr(1,1 + powerMatInit(3,ii)*gammaRRT(ii) + powerMatInit(1,ii)*gammaSRT(ii))...
                          - rel_entr(1,1 + powerMatInit(1,ii)*gammaSDT(ii) + powerMatInit(2,ii)*gammaSDT(ii))...
                          + gammaRRT(ii)/(1 + powerMatInit(3,ii)*gammaRRT(ii) + powerMatInit(1,ii)*gammaSRT(ii))*(powerMat(3,ii)-powerMatInit(3,ii))...
                          + gammaSRT(ii)/(1 + powerMatInit(3,ii)*gammaRRT(ii) + powerMatInit(1,ii)*gammaSRT(ii))*(powerMat(1,ii)-powerMatInit(1,ii))...
                          + gammaSDT(ii)/(1 + powerMatInit(1,ii)*gammaSDT(ii) + powerMatInit(2,ii)*gammaSDT(ii))*(powerMat(1,ii)-powerMatInit(1,ii))...
                          + gammaSDT(ii)/(1 + powerMatInit(1,ii)*gammaSDT(ii) + powerMatInit(2,ii)*gammaSDT(ii))*(powerMat(2,ii)-powerMatInit(2,ii))...
                          + cons(ii);
                   end
                   minimize obj(2*N+1)
                   subject to
                       cons(2*N+1) <= 0
                       sum(sum(powerMat(1:2,:))) <= 2*P_S
                       sum(sum(powerMat(3,:))) <= 2*P_R
                       powerMat >= 0
               cvx_end
               powerMatInit = powerMat;
               optvalEor = cvx_optval - cvx_optvalInit;
               cvx_optvalInit = cvx_optval;
        end
%         [powerMat3,cvx_optval3] = scheme1_SCP(powerMatInit3,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);       
%         [powerMat4,cvx_optval4] = scheme1_SCP(powerMatInit4,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);       
%         cvx_optval = min([cvx_optval2,cvx_optval3,cvx_optval4]);
        rate(kk) = -cvx_optval2/log(2)/N/2;
    end
    rateAve = rateAve + rate;
    toc
end
rateAve = rateAve./T;













































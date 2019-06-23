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
T       = 1; % channel realization number

rate = 0;
rateAve = zeros(length(Rho),1);

for kk = 13
    tt = 0;
    while(tt<T)
        kk
        tt
        Hsr     = omega_sr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
        Hrr     = omega_rr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
        Hrd     = omega_rd/sqrt(2)*(randn(N,1)+1i*randn(N,1));
        Hsd     = omega_sd/sqrt(2)*(randn(N,1)+1i*randn(N,1));    
        gammaSRT = abs(Hsr).^2;
        gammaRDT = abs(Hrd).^2;
        gammaSDT = abs(Hsd).^2;
        gammaRRT = Rho(kk)*abs(Hrr).^2;
        
        gammaSRT = [1.313755645961485;1.890351625082766;0.980904518015875;0.626475536012609];
        gammaRDT = [0.113738932920653;0.324317426297855;0.216468700589692;1.253387322360803];
        gammaSDT = [0.006706091425047;0.007882029639352;0.003954711374731;0.012078472643101];
        gammaRRT = [6.956275097126855;4.843780079986045;6.608475818948081;3.845590925544955];
        
        % initialized points
        powerMatInit2 = powerMatInitTwo(P_S,P_R,N,gammaSRT,gammaRDT);
        powerMatInit3 = powerMatInitThree(P_S,N,gammaSDT);
        powerMatInit4 = powerMatInitFour(P_S,P_R,N,gammaSRT,gammaRDT);
        
        timeRatioInit = 0.2;
        powerMatInit = powerMatInit2;
        
        timeRatio = timeRatioInit;
        powerMat = powerMatInit;
        
        % check feasibility
%         feasId2 = scheme2_SCP_CheckFeasibility(timeRatioInit,powerMatInit2,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
%         feasId3 = scheme2_SCP_CheckFeasibility(powerMatInit3,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
%         feasId4 = scheme2_SCP_CheckFeasibility(powerMatInit4,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
        cvx_begin quiet
        cvx_solver Mosek_2
           variable powerMat(3,2*N)
           expression cons
           cons = 0;
           for ii = 1:N
              cons =  timeRatio.*(-log(1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
                  - log(1 + powerMat(3,2*ii-1)*gammaRDT(ii) + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
                  + log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))...
                  + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))...
                  + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(3,2*ii-1)-powerMatInit(3,2*ii-1))...
                  + gammaSRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
                  + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
                  + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(2,2*ii-1)-powerMatInit(2,2*ii-1)))...
                  +(1-timeRatio)*(-log(1 + powerMat(3,2*ii)*gammaRRT(ii))...
                  - log(1 + powerMat(3,2*ii)*gammaRDT(ii) + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
                  + log(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))...
                  + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))...
                  + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(3,2*ii)-powerMatInit(3,2*ii))...
                  + gammaSRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
                  + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
                  + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(2,2*ii)-powerMatInit(2,2*ii)))...
                  + cons;
           end
           minimize 1
           subject to
               cons <= 0
               sum(sum(timeRatio.*powerMat(1:2,1:2:2*N) + (1-timeRatio).*powerMat(1:2,2:2:2*N))) <= P_S
               sum(sum(timeRatio.*powerMat(3,1:2:2*N) + (1-timeRatio).*powerMat(3,2:2:2*N))) <= P_R
               powerMat >= 0
        cvx_end
        if cvx_status(1) == 'S'
            feasId2 = 1;
        else
            feasId2 = 0;
        end
        if feasId2 == 0% || feasId3 == 0 || feasId4 == 0
            continue; % if one is not feasible, start a new channel realization
        end
        
        % rate optimization
        rateInit = 0;
        rateEpi = 10;
        iter = 0;
        iterMax = 30;
        tic
%         while(rateEpi > 1e-3)
            cvx_begin quiet
            cvx_solver Mosek_2
               variable powerMat(3,2*N)
               expression obj
               expression cons
               obj = 0;
               cons = 0;
               for ii = 1:N
                  obj = timeRatio*(log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))...
                      + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))...
                      + rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRRT(ii) + powerMat(1,2*ii-1)*gammaSRT(ii))...
                      + rel_entr(1,1 + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
                      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))*(powerMat(3,2*ii-1) - powerMatInit(3,2*ii-1))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1) - powerMatInit(1,2*ii-1)))...
                      + (1-timeRatio)*(log(1 + powerMatInit(3,2*ii)*gammaRRT(ii))...
                      + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii))...
                      + rel_entr(1,1 + powerMat(3,2*ii)*gammaRRT(ii) + powerMat(1,2*ii)*gammaSRT(ii))...
                      + rel_entr(1,1 + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
                      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii))*(powerMat(3,2*ii) - powerMatInit(3,2*ii))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii) - powerMatInit(1,2*ii)))...
                      + obj;
                  cons =  timeRatio*(rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
                      + rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRDT(ii) + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
                      + log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))...
                      + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))...
                      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(3,2*ii-1)-powerMatInit(3,2*ii-1))...
                      + gammaSRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(2,2*ii-1)-powerMatInit(2,2*ii-1)))...
                      +(1-timeRatio)*(rel_entr(1,1 + powerMat(3,2*ii)*gammaRRT(ii))...
                      + rel_entr(1,1 + powerMat(3,2*ii)*gammaRDT(ii) + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
                      + log(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))...
                      + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))...
                      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(3,2*ii)-powerMatInit(3,2*ii))...
                      + gammaSRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
                      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(2,2*ii)-powerMatInit(2,2*ii)))...
                      + cons;
%                   obj = timeRatio*(log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))...
%                       + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))...
%                       -log(1 + powerMat(3,2*ii-1)*gammaRRT(ii) + powerMat(1,2*ii-1)*gammaSRT(ii))...
%                       -log(1 + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
%                       + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))*(powerMat(3,2*ii-1) - powerMatInit(3,2*ii-1))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1) - powerMatInit(1,2*ii-1)))...
%                       + (1-timeRatio)*(log(1 + powerMatInit(3,2*ii)*gammaRRT(ii))...
%                       + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii))...
%                       -log(1 + powerMat(3,2*ii)*gammaRRT(ii) + powerMat(1,2*ii)*gammaSRT(ii))...
%                       -log(1 + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
%                       + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii))*(powerMat(3,2*ii) - powerMatInit(3,2*ii))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii) - powerMatInit(1,2*ii)))...
%                       + obj;
%                   cons =  timeRatio*(-log(1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
%                       -log(1 + powerMat(3,2*ii-1)*gammaRDT(ii) + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
%                       + log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))...
%                       + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))...
%                       + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(3,2*ii-1)-powerMatInit(3,2*ii-1))...
%                       + gammaSRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii) + powerMatInit(2,2*ii-1)*gammaSDT(ii))*(powerMat(2,2*ii-1)-powerMatInit(2,2*ii-1)))...
%                       +(1-timeRatio)*(-log(1 + powerMat(3,2*ii)*gammaRRT(ii))...
%                       -log(1 + powerMat(3,2*ii)*gammaRDT(ii) + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
%                       + log(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))...
%                       + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))...
%                       + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(3,2*ii)-powerMatInit(3,2*ii))...
%                       + gammaSRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
%                       + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii) + powerMatInit(2,2*ii)*gammaSDT(ii))*(powerMat(2,2*ii)-powerMatInit(2,2*ii)))...
%                       + cons;
               end
               minimize obj
               subject to
                   cons <= 0
                   sum(sum(timeRatio.*powerMat(1:2,1:2:2*N) + (1-timeRatio).*powerMat(1:2,2:2:2*N))) <= P_S
                   sum(sum(timeRatio.*powerMat(3,1:2:2*N) + (1-timeRatio).*powerMat(3,2:2:2*N))) <= P_R
                   powerMat >= 0
            cvx_end
%             if isnan(cvx_optval) == 1
%                 a = 1;
%             end
%         end
       
        
%         tic
%         [~,rate3,iter3] = scheme2_SCP(powerMatInit3,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);       
%         toc
%         
%         tic 
%         [~,rate4,iter4] = scheme2_SCP(powerMatInit4,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);       
%         toc
        
        tt = tt + 1;
    end
    toc
end

    












































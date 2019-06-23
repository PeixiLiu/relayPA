function [rate,timeRatio,powerMat,iterCount] = alterOptim3(timeRatioInit, powerMatInit,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT)
%alterOptim Alternately optimization algorithm
%   
rateInit = 0;
rateEpi = 10;
iter = 1;
iterMax = 100;
while(rateEpi > 1e-3 && iter <= iterMax)
    powerMat = scheme3_SCP(timeRatioInit, powerMatInit,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
    timeRatio = scheme3_timeRatio(powerMat,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
    rate = scheme3_rate(timeRatio,powerMat,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT);
    rateEpi = abs(rate - rateInit);
    rateInit = rate;
    timeRatioInit = timeRatio;
    powerMatInit = powerMat;
    iter = iter + 1;
end
iterCount = iter;
end


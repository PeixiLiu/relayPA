function [feasibilityState] = scheme3_SCP_CheckFeasibility(timeRatio, powerMatInit,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT)
%scheme3_SCP SCP algorithm aa
%   powerMat: optimal power allocation
%   cvx_optval: optimal value of cvx

cvx_begin quiet
cvx_solver 
   variable powerMat(3,2*N)
   expression cons
   cons = 0;
   for ii = 1:N
      cons =  timeRatio*(rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
          + rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRDT(ii) + powerMat(1,2*ii-1)*gammaSDT(ii))...
          + log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))...
          + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))...
          + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(3,2*ii-1)-powerMatInit(3,2*ii-1))...
          + gammaSRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii) + powerMatInit(1,2*ii-1)*gammaSRT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1))...
          + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1)-powerMatInit(1,2*ii-1)))...
          +(1-timeRatio)*(rel_entr(1,1 + powerMat(3,2*ii)*gammaRRT(ii))...
          + rel_entr(1,1 + powerMat(3,2*ii)*gammaRDT(ii) + powerMat(1,2*ii)*gammaSDT(ii))...
          + log(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))...
          + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii))...
          + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(3,2*ii)-powerMatInit(3,2*ii))...
          + gammaSRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii) + powerMatInit(1,2*ii)*gammaSRT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii))...
          + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii)-powerMatInit(1,2*ii)))...
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
    feasibilityState = 1;
else
    feasibilityState = 0;
end
end


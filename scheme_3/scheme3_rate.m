function [rate] = scheme3_rate(timeRatio,powerMat,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT)
%scheme2_SCP SCP algorithm aa
%   powerMat: optimal power allocation
%   cvx_optval: optimal value of cvx
obj = 0;
for ii = 1:N
  obj = timeRatio*(log(1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
      + log(1 + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(3,2*ii-1)*gammaRDT(ii))...
      + rel_entr(1,1 + powerMat(3,2*ii-1)*gammaRRT(ii) + powerMat(1,2*ii-1)*gammaSRT(ii))...
      + rel_entr(1,1 + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii) + powerMat(3,2*ii-1)*gammaRDT(ii)))...
      + (1-timeRatio)*(log(1 + powerMat(3,2*ii)*gammaRRT(ii))...
      + log(1 + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(3,2*ii)*gammaRDT(ii))...
      + rel_entr(1,1 + powerMat(3,2*ii)*gammaRRT(ii) + powerMat(1,2*ii)*gammaSRT(ii))...
      + rel_entr(1,1 + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii) + powerMat(3,2*ii)*gammaRDT(ii)))...
      + obj;
end
rate = -obj./N./log(2);
end



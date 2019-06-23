function [powerMat] = scheme2_SCP(timeRatio, powerMatInit,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT)
%scheme2_SCP SCP algorithm aa
%   powerMat: optimal power allocation
%   cvx_optval: optimal value of cvx

%% YALMIP
% Define variables
powerMat = sdpvar(3,2*N);

% Define constraints
cons = 0;
for ii = 1:N
  cons =  timeRatio*(-log(1 + powerMat(3,2*ii-1)*gammaRRT(ii))...
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
Constraints = [cons <= 0,...
    sum(sum(timeRatio.*powerMat(1:2,1:2:2*N) + (1-timeRatio).*powerMat(1:2,2:2:2*N))) <= P_S,...
    sum(sum(timeRatio.*powerMat(3,1:2:2*N) + (1-timeRatio).*powerMat(3,2:2:2*N))) <= P_R,...
    powerMat >= 0];

% Define an objective
obj = 0;
for ii = 1:N
  obj = timeRatio*(log(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))...
      + log(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))...
      - log(1 + powerMat(3,2*ii-1)*gammaRRT(ii) + powerMat(1,2*ii-1)*gammaSRT(ii))...
      - log(1 + powerMat(1,2*ii-1)*gammaSDT(ii) + powerMat(2,2*ii-1)*gammaSDT(ii))...
      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii-1)*gammaRRT(ii))*(powerMat(3,2*ii-1) - powerMatInit(3,2*ii-1))...
      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii-1)*gammaSDT(ii))*(powerMat(1,2*ii-1) - powerMatInit(1,2*ii-1)))...
      + (1-timeRatio)*(log(1 + powerMatInit(3,2*ii)*gammaRRT(ii))...
      + log(1 + powerMatInit(1,2*ii)*gammaSDT(ii))...
      - log(1 + powerMat(3,2*ii)*gammaRRT(ii) + powerMat(1,2*ii)*gammaSRT(ii))...
      - log(1 + powerMat(1,2*ii)*gammaSDT(ii) + powerMat(2,2*ii)*gammaSDT(ii))...
      + gammaRRT(ii)/(1 + powerMatInit(3,2*ii)*gammaRRT(ii))*(powerMat(3,2*ii) - powerMatInit(3,2*ii))...
      + gammaSDT(ii)/(1 + powerMatInit(1,2*ii)*gammaSDT(ii))*(powerMat(1,2*ii) - powerMatInit(1,2*ii)))...
      + obj;
end
Objective = obj;

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','mosek');

% Solve the problem
sol = optimize(Constraints,Objective,options);
cvx_optval = Objective;
powerMat = value(powerMat);
% Analyze error flags
if sol.problem ~= 0
    disp('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end

end

function [timeRatio] = scheme1_timeRatio(powerMat,P_S,P_R,N,gammaSRT,gammaRDT,gammaSDT,gammaRRT)
%scheme1_SCP timeRatio

%%
% Define variables
timeRatio = sdpvar(1,1);

% Define constraints
cons = 0;
for ii = 1:N
  cons =  timeRatio*(-log(1 + powerMat(2,2*ii-1)*gammaSRT(ii) + powerMat(3,2*ii-1)*gammaRRT(ii))...
      - log(1+ (powerMat(1,2*ii-1) + powerMat(2,2*ii-1))*gammaSDT(ii) + powerMat(3,2*ii-1)*gammaRDT(ii))...
      + log(1 + (powerMat(1,2*ii-1) + powerMat(2,2*ii-1))*gammaSRT(ii) + powerMat(3,2*ii-1)*gammaRRT(ii))...
      + log(1 + (powerMat(1,2*ii-1) + powerMat(2,2*ii-1))*gammaSDT(ii)))...
      +(1-timeRatio)*(-log(1 + powerMat(2,2*ii)*gammaSRT(ii) + powerMat(3,2*ii)*gammaRRT(ii))...
      - log(1+ (powerMat(1,2*ii) + powerMat(2,2*ii))*gammaSDT(ii) + powerMat(3,2*ii)*gammaRDT(ii))...
      + log(1 + (powerMat(1,2*ii) + powerMat(2,2*ii))*gammaSRT(ii) + powerMat(3,2*ii)*gammaRRT(ii))...
      + log(1 + (powerMat(1,2*ii) + powerMat(2,2*ii))*gammaSDT(ii)))...
      + cons;
end
Constraints = [cons <= 0,...
    sum(sum(timeRatio.*powerMat(1:2,1:2:2*N) + (1-timeRatio).*powerMat(1:2,2:2:2*N))) <= P_S,...
    sum(sum(timeRatio.*powerMat(3,1:2:2*N) + (1-timeRatio).*powerMat(3,2:2:2*N))) <= P_R,...
    0 <= timeRatio <=1];

% Define an objective
obj = 0;
for ii = 1:N
  obj = timeRatio.*(-log(1 + powerMat(1,2*ii-1)*gammaSRT(ii) + powerMat(2,2*ii-1)*gammaSRT(ii)...
      + powerMat(3,2*ii-1)*gammaRRT(ii)) - log(1 + powerMat(2,2*ii-1)*gammaSDT(ii))...
      + log(1 + powerMat(2,2*ii-1)*gammaSRT(ii) + powerMat(3,2*ii-1)*gammaRRT(ii)))...
      + (1 - timeRatio).*(-log(1 + powerMat(1,2*ii)*gammaSRT(ii) + powerMat(2,2*ii)*gammaSRT(ii)...
      + powerMat(3,2*ii)*gammaRRT(ii)) - log(1 + powerMat(2,2*ii)*gammaSDT(ii))...
      + log(1 + powerMat(2,2*ii)*gammaSRT(ii) + powerMat(3,2*ii)*gammaRRT(ii)))...
      + obj;
end
Objective = obj;

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','mosek');

% Solve the problem
sol = optimize(Constraints,Objective,options);
timeRatio = value(timeRatio);
% Analyze error flags
if sol.problem ~= 0
    disp('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end
if timeRatio > 1
    timeRatio = 2 - timeRatio; % in case value a little bigge than 1 due to precision of matlab
end
if timeRatio < 0
   timeRatio = abs(timeRatio); % in case small negative value due to precision of matlab 
end
end

function [powerMatInit] = powerMatInitThree(P_S,N,gammaSD)
%powerMatInitOne power initialization method three (waterfilling of direct link)
% P_S  : Source power constraint
% P_R  : Relay power constraint
% N    : number of subcarrier

powerSDWF = walter_filling(gammaSD,N,P_S);

powerMatInit(1,:)=0;
powerMatInit(2,1:2:2*N-1) = powerSDWF; % pSR, pSD, pR
powerMatInit(2,2:2:2*N) = powerSDWF;
powerMatInit(3,:) = 0;

end


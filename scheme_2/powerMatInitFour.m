function [powerMatInit] = powerMatInitFour(P_S,P_R,N,gammaSR,gammaRD)
%powerMatInitOne power initialization method four (waterfilling of relay link)
% P_S  : Source power constraint
% P_R  : Relay power constraint
% N    : number of subcarrier

powerSRWF = walter_filling(gammaSR,N,P_S);
powerRDWF = walter_filling(gammaRD,N,P_R);

powerMatInit(1,1:2:2*N-1) = powerSRWF; % pSR, pSD, pR
powerMatInit(1,2:2:2*N) = 0;
powerMatInit(2,:)=0;
powerMatInit(3,1:2:2*N-1) = 0;
powerMatInit(3,2:2:2*N) = powerRDWF;

end


function [powerMatInit] = powerMatInitTwo(P_S,P_R,N,gammaSR,gammaRD,gammaRR,gammaSD)
%powerMatInitOne power initialization method two (waterfilling of relay link)
% P_S  : Source power constraint
% P_R  : Relay power constraint
% N    : number of subcarrier

% powerSRWF = walter_filling(gammaSR,N,P_S);
% powerRDWF = walter_filling(gammaRD,N,P_R);
% 
% powerMatInit(1,1:2:2*N-1) = powerSRWF; % pSR, pSD, pR
% powerMatInit(1,2:2:2*N) = powerSRWF;
% powerMatInit(2,:)=0;
% powerMatInit(3,1:2:2*N-1) = powerRDWF;
% powerMatInit(3,2:2:2*N) = powerRDWF;


% Hsr     = omega_sr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
% Hrr     = omega_rr/sqrt(2)*(randn(N,1)+1i*randn(N,1));
% Hrd     = omega_rd/sqrt(2)*(randn(N,1)+1i*randn(N,1));
% Hsd     = omega_sd/sqrt(2)*(randn(N,1)+1i*randn(N,1));    
A       = gammaSR;
C       = gammaRD;
D       = gammaSD;
B       = gammaRR;

%half-duplex

X1 = walter_filling(A,N,P_S);
Y1 = walter_filling(C,N,P_R);


% SCP iteration parameters
Num     = 10;
Epsilon = 1e-2;
num     = 0;
epsilon = 2*Epsilon;

X       = X1;
Y       = Y1;
while( num<=Num && epsilon>Epsilon )
    x = X;
    y = Y;
    X = Source_Allocation(A,B,C,D,N,P_S,X,Y);
    Y = Source_Allocation(C,D,A,B,N,P_R,Y,X);
    %Y = Relay_Allocation(A,B,C,D,N,P_R(i),X,Y);
    epsilon = norm([x-X;y-Y])/(2*N*norm([X;Y]));
    num = num + 1;
end
powerMatInit(1,1:2:2*N-1) = X'; % pSR, pSD, pR
powerMatInit(1,2:2:2*N) = X';
powerMatInit(2,:)=0;
powerMatInit(3,1:2:2*N-1) = Y';
powerMatInit(3,2:2:2*N) = Y';
%     X_A(k,:) = X;
%     Y_A(k,:) = Y;
%     R_temp1 = Capacity(A,B,N,X,Y);
%     R_temp2 = Capacity(C,D,N,Y,X);
end

function X = Source_Allocation(A,B,C,D,N,P_S,X,Y)
Num_S       = 20;
Epsilon_S   = 1e-3;
num_S       = 0;
epsilon_S   = 2*Epsilon_S;
while (num_S<=Num_S && epsilon_S>Epsilon_S)
    x = X;
    X = Source_2D_Search(A,B,C,D,N,P_S,X,Y);
    epsilon_S = norm(X-x)/norm(X)/N;
    num_S     = num_S + 1;
end
end

function X_new = Source_2D_Search(A,B,C,D,N,P_S,X_old,Y)
e          = exp(1);
threshold  = 1/N*sum(log2(1+A.*P_S./N./(1+B.*Y)));
phi        = log2(e)/N*(A./(1+A.*X_old+B.*Y)-D./(1+D.*X_old+C.*Y)+D./(1+D.*X_old));
%==========================================================================
% Case I: \lambda ~= 0, \upsilon = 0
lambda_temp1 = log2(e)/N*(A./(1+A.*P_S+B.*Y))./phi;
lambda_temp2 = log2(e)/N*(A./(1+B.*Y))./phi;
lambda_min   = max(lambda_temp1);
lambda_max   = max(lambda_temp2);
X = Solve_X(A,B,Y,N,lambda_min,0,phi);
hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi);
if hk_X>0 % The first constraint can be achieved
    X       = zeros(N,1);
    Epsilon1 = 1e-2*threshold;
    mu      = 2*Epsilon1;
    while(mu>Epsilon1) 
        lambda  = (lambda_min + lambda_max)/2;
        for n = 1:N
            % calculate X for a given lambda
            if lambda >= lambda_temp2(n)
                X(n) = 0;
            else
                X(n) = Solve_X(A(n),B(n),Y(n),N,lambda,0,phi(n));
            end
        end
        % bisection search
        hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi);
        mu    = abs(hk_X);
        if hk_X > 0
            lambda_min = lambda;
        else
            lambda_max = lambda;
        end 
    end
    if sum(X)<P_S % The second constraint can be achieved
        X_new = X;
        return;
    end
end
%==========================================================================
% Case II: \lambda = 0, \upsilon ~= 0
upsilon_temp1 = log2(e)/N*(A./(1+A.*P_S+B.*Y));
upsilon_temp2 = log2(e)/N*(A./(1+B.*Y));
upsilon_min   = max(upsilon_temp1);
upsilon_max   = max(upsilon_temp2);
X = Solve_X(A,B,Y,N,0,upsilon_min,phi);
if sum(X)>P_S % The second constraint can be achieved
    X       = zeros(N,1);
    Epsilon2 = 1e-5*P_S;
    mu      = 2*Epsilon2;
    while(abs(mu-P_S)>Epsilon2) 
        upsilon  = (upsilon_min + upsilon_max)/2;
        for n = 1:N
            % calculate X for a given lambda
            if upsilon >= upsilon_temp2(n)
                X(n) = 0;
            else
                X(n) = Solve_X(A(n),B(n),Y(n),N,0,upsilon,0);
            end
        end
        % bisection search
        mu    = sum(X);
        if mu > P_S
            upsilon_min = upsilon;
        else
            upsilon_max = upsilon;
        end 
    end
    hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi);
    if hk_X<0 % The first constraint can be achieved
        X_new = X;
        return;
    end
else
    upsilon = 0;
end
%==========================================================================
% Case III: \lambda ~= 0, \upsilon ~= 0
It_Num  = 200;
It_Num1 = 200;
It_Num2 = 200;
it      = 0;
while(1)
    it = it + 1;
    if it == It_Num
        break;
    end
    %subcase I: search \lambda ============================================
    lambda_temp1 = (log2(e)/N*(A./(1+A.*P_S+B.*Y))-upsilon)./phi;
    lambda_temp2 = (log2(e)/N*(A./(1+B.*Y))-upsilon)./phi;
    lambda_min   = max(lambda_temp1);
    lambda_max   = max(lambda_temp2);
    X       = zeros(N,1);
    Epsilon1 = 1e-2*threshold;
    mu      = 2*Epsilon1;
    it1     = 0;
    while(mu>Epsilon1) 
        it1 = it1 + 1;
        if it1 == It_Num1
            break;
        end
        lambda  = (lambda_min + lambda_max)/2;
        for n = 1:N
            % calculate X for a given lambda
            if lambda >= lambda_temp2(n)
                X(n) = 0;
            else
                X(n) = Solve_X(A(n),B(n),Y(n),N,lambda,upsilon,phi(n));
            end
        end
        % bisection search
        hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi);
        mu    = abs(hk_X);
        if hk_X > 0
            lambda_min = lambda;
        else
            lambda_max = lambda;
        end 
    end
    %subcase II: search \upsilon ==========================================
    upsilon_temp1 = log2(e)/N*(A./(1+A.*P_S+B.*Y))-lambda*phi;
    upsilon_temp2 = log2(e)/N*(A./(1+B.*Y))-lambda*phi;
    upsilon_min   = max(upsilon_temp1);
    upsilon_max   = max(upsilon_temp2);
    X       = zeros(N,1);
    Epsilon2 = 1e-5*P_S;
    mu      = 2*Epsilon2;
    it2     = 0;
    while(abs(mu-P_S)>Epsilon2)
        it2 = it2 + 1;
        if it2== It_Num2
            break;
        end
        upsilon  = (upsilon_min + upsilon_max)/2;
        for n = 1:N
            % calculate X for a given lambda
            if upsilon >= upsilon_temp2(n)
                X(n) = 0;
            else
                X(n) = Solve_X(A(n),B(n),Y(n),N,lambda,upsilon,phi(n));
            end
        end
        % bisection search
        mu    = sum(X);
        if mu > P_S
            upsilon_min = upsilon;
        else
            upsilon_max = upsilon;
        end
    end
    hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi);
    if abs(hk_X)<Epsilon1 % The first constraint can be achieved
        X_new = X;
        return;
    end
end
X_new  = X;
end

function X = Solve_X(A,B,Y,N,lambda,upsilon,phi)
e = exp(1);
X = (log2(e)/N*A./(lambda*phi+upsilon)-1-B.*Y)./A;
X = X.*(X>0);
end

function hk_X = HK_X(A,B,C,D,X,X_old,Y,N,phi)
% temp1 = log2(1+A.*X_old./(1+B.*Y))-log2(1+C.*Y./(1+D.*X_old));
% temp2 = phi.*(X-X_old);
% hk_X = sum(temp1)/N+sum(temp2);
temp1 = log2(1+A.*X./(1+B.*Y))-log2(1+C.*Y./(1+D.*X));
hk_X = sum(temp1)/N;
end


function R = Capacity(A,B,N,X,Y)
temp1  = log2(1+A.*X./(1+B.*Y));
R     = 1/N*sum(temp1);
end

%==========================================================================
%half-duplex 
function X = walter_filling(A,N,P)
e            = exp(1);
lambda_temp1 = log2(e)/N*(A./(1+A.*P));
lambda_temp2 = log2(e)/N*A;
lambda_min   = max(lambda_temp1);
lambda_max   = max(lambda_temp2);
X            = zeros(N,1);
Epsilon      = 1e-5*P;
mu           = 2*Epsilon;
while(mu>Epsilon)
    lambda  = (lambda_min + lambda_max)/2;
    for n = 1:N
        % calculate X for a given lambda
        if lambda >= lambda_temp2(n)
            X(n) = 0;
        else
            X(n) = log2(e)/N/lambda-1/A(n);
        end
    end
    % bisection search
    mu    = abs(sum(X)-P);
    if sum(X) > P
        lambda_min = lambda;
    else
        lambda_max = lambda;
    end
end
end

function R = Capacity_HD(A,X,N)
temp1  = log2(1+A.*X);
R1     = 1/N*sum(temp1);
R      = R1/2;
end
%==========================================================================
%OPtimal_Power_Allocation with symbol wise processing
% full-duplex, with direct link
function [X,Y] = FD_PA_Individual_Constraint(A,B,C,D,P_S,P_R,N)
X         = zeros(N,1);
Y         = zeros(N,1);
temp_phi1 = C.*D-A.*B;
temp_phi2 = 1:N;
phi1      = temp_phi2(temp_phi1>=0);
phi2      = temp_phi2(temp_phi1<0);
y_temp    = fn(A,B,C,D,P_S);
x_temp    = gn(A,B,C,D,P_R);
% lambda ~= 0, upsilon = 0
upsilon    = 0;
Lambda_Max = zeros(N,1);
Lambda_Min = zeros(N,1);
for n = phi1
    Lambda_Max(n) = Fn_lambda(A(n),B(n),C(n),D(n),0,N,upsilon);
    Lambda_Min(n) = Fn_lambda(A(n),B(n),C(n),D(n),P_S,N,upsilon);
end
for n = phi2
    Lambda_Max(n) = Gn_lambda(A(n),B(n),C(n),D(n),0,N,upsilon);
    Lambda_Min(n) = Gn_lambda(A(n),B(n),C(n),D(n),y_temp(n),N,upsilon);
end
lambda_max = max(Lambda_Max);
lambda_min = max(max(Lambda_Min),0);
Epsilon1    = 1e-5*P_S;
epsilon1    = 2*Epsilon1;
it1        = 0;
It_Num1    = 100;
while (abs(epsilon1)>Epsilon1)
    lambda = (lambda_max + lambda_min)/2;
    for n = phi1
        if lambda >  Lambda_Max(n)
            X(n) = 0;
        else
            X(n) = Fn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,P_S);
        end
        Y(n) = fn(A(n),B(n),C(n),D(n),X(n));
    end
    for n = phi2
        if lambda >  Lambda_Max(n)
            Y(n) = 0;
        else
            Y(n) = Gn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,y_temp(n));
        end
        X(n) = gn(A(n),B(n),C(n),D(n),Y(n));
    end
    epsilon1 = sum(X) - P_S;
    if epsilon1 > 0
        lambda_min = lambda;
    else
        lambda_max = lambda;
    end
    it1 = it1 + 1;
    if it1 == It_Num1
        break;
    end
end
if sum(Y) < P_R
    return;
end
% lambda = 0, upsilon ~= 0
lambda    = 0;
Upsilon_Max = zeros(N,1);
Upsilon_Min = zeros(N,1);
for n = phi1
    Upsilon_Max(n) = Fn_upsilon(A(n),B(n),C(n),D(n),0,N,lambda);
    Upsilon_Min(n) = Fn_upsilon(A(n),B(n),C(n),D(n),x_temp(n),N,lambda);
end
for n = phi2
    Upsilon_Max(n) = Gn_upsilon(A(n),B(n),C(n),D(n),0,N,lambda);
    Upsilon_Min(n) = Gn_upsilon(A(n),B(n),C(n),D(n),P_R,N,lambda);
end
upsilon_max = max(Upsilon_Max);
upsilon_min = max(max(Upsilon_Min),0);
Epsilon2    = 1e-5*P_R;
epsilon2    = 2*Epsilon2;
it2        = 0;
It_Num2    = 100;
while (abs(epsilon2)>Epsilon2)
    upsilon = (upsilon_max + upsilon_min)/2;
    for n = phi1
        if upsilon > Upsilon_Max(n)
            X(n) = 0;
        else
            X(n) = Fn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,x_temp(n));
        end
        Y(n) = fn(A(n),B(n),C(n),D(n),X(n));
    end
    for n = phi2
        if upsilon > Upsilon_Max(n)
            Y(n) = 0;
        else
            Y(n) = Gn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,P_R);
        end
        X(n) = gn(A(n),B(n),C(n),D(n),Y(n));
    end
    epsilon2 = sum(Y) - P_R;
    if epsilon2 > 0
        upsilon_min = upsilon;
    else
        upsilon_max = upsilon;
    end
    it2 = it2 + 1;
    if it2 == It_Num2
        break;
    end
end
if sum(X) < P_S
    return;
end
% lambda ~= 0, upsilon ~= 0
Epsilon3  = 1e-4*P_S;
It_Num    = 200;
it        = 0;
while (1)
    it = it + 1;
    if it == It_Num
        break;
    end
	% search lambda when upsilon is fixed
    Lambda_Max = zeros(N,1);
    Lambda_Min = zeros(N,1);
    for n = phi1
        Lambda_Max(n) = Fn_lambda(A(n),B(n),C(n),D(n),0,N,upsilon);
        Lambda_Min(n) = Fn_lambda(A(n),B(n),C(n),D(n),P_S,N,upsilon);
    end
    for n = phi2
        Lambda_Max(n) = Gn_lambda(A(n),B(n),C(n),D(n),0,N,upsilon);
        Lambda_Min(n) = Gn_lambda(A(n),B(n),C(n),D(n),y_temp(n),N,upsilon);
    end
    lambda_max = max(Lambda_Max);
    lambda_min = max(max(Lambda_Min),0);
    Epsilon1   = 1e-5*P_S;
    epsilon1   = 2*Epsilon1;
    it1        = 0;
    It_Num1    = 100;
    while (abs(epsilon1)>Epsilon1)
        lambda = (lambda_max + lambda_min)/2;
        for n = phi1
            if lambda >  Lambda_Max(n)
                X(n) = 0;
            else
                X(n) = Fn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,P_S);
            end
            Y(n) = fn(A(n),B(n),C(n),D(n),X(n));
        end
        for n = phi2
            if lambda >  Lambda_Max(n)
                Y(n) = 0;
            else
                Y(n) = Gn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,y_temp(n));
            end
            X(n) = gn(A(n),B(n),C(n),D(n),Y(n));
        end
        epsilon1 = sum(X) - P_S;
        if epsilon1 > 0
            lambda_min = lambda;
        else
            lambda_max = lambda;
        end
        it1 = it1 + 1;
        if it1 == It_Num1
            break;
        end
    end
    % search upsilon when lambda is fixed
    Upsilon_Max = zeros(N,1);
    Upsilon_Min = zeros(N,1);
    for n = phi1
        Upsilon_Max(n) = Fn_upsilon(A(n),B(n),C(n),D(n),0,N,lambda);
        Upsilon_Min(n) = Fn_upsilon(A(n),B(n),C(n),D(n),x_temp(n),N,lambda);
    end
    for n = phi2
        Upsilon_Max(n) = Gn_upsilon(A(n),B(n),C(n),D(n),0,N,lambda);
        Upsilon_Min(n) = Gn_upsilon(A(n),B(n),C(n),D(n),P_R,N,lambda);
    end
    upsilon_max = max(Upsilon_Max);
    upsilon_min = max(max(Upsilon_Min),0);
    Epsilon2    = 1e-5*P_R;
    epsilon2    = 2*Epsilon2;
    it2        = 0;
    It_Num2    = 100;
    while (abs(epsilon2)>Epsilon2)
        upsilon = (upsilon_max + upsilon_min)/2;
        for n = phi1
            if upsilon > Upsilon_Max(n)
                X(n) = 0;
            else
                X(n) = Fn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,x_temp(n));
            end
            Y(n) = fn(A(n),B(n),C(n),D(n),X(n));
        end
        for n = phi2
            if upsilon > Upsilon_Max(n)
                Y(n) = 0;
            else
                Y(n) = Gn_inv(A(n),B(n),C(n),D(n),N,lambda,upsilon,P_R);
            end
            X(n) = gn(A(n),B(n),C(n),D(n),Y(n));
        end
        epsilon2 = sum(Y) - P_R;
        if epsilon2 > 0
            upsilon_min = upsilon;
        else
            upsilon_max = upsilon;
        end
        it2 = it2 + 1;
        if it2 == It_Num2
            break;
        end
    end
    if abs(sum(X)-P_S) < Epsilon3
        return;
    end
end
end

function y = Fn(A,B,C,D,x,N)
temp    = C.^2 + 4*A.*B.*C*x + 4*A.*B.*C.*D*x^2;
Delta   = sqrt(temp);
f       = (-C+Delta)./(2*B.*C);
f1      = (A+2*A.*D*x)./Delta;
y       = (A.*(1+B.*f)-A.*B.*f1*x)./((1+B.*f).^2+A.*(1+B.*f)*x)*log2(exp(1))/N;
%y       = lambda + upsilon*f1 - y1;
end

function y = Fn_lambda_upsilon(A,B,C,D,x,lambda,upsilon)
temp    = C.^2 + 4*A.*B.*C*x + 4*A.*B.*C.*D*x^2;
Delta   = sqrt(temp);
f1      = (A+2*A.*D*x)./Delta;
y       = lambda + upsilon*f1;
end

function y = Fn_lambda(A,B,C,D,x,N,upsilon)
temp    = C.^2 + 4*A.*B.*C*x + 4*A.*B.*C.*D*x^2;
Delta   = sqrt(temp);
f       = (-C+Delta)./(2*B.*C);
f1      = (A+2*A.*D*x)./Delta;
y1      = (A.*(1+B.*f)-A.*B.*f1*x)./((1+B.*f).^2+A.*(1+B.*f)*x)*log2(exp(1))/N;
y       = y1 - upsilon*f1;
end

function y = fn(A,B,C,D,x)
temp    = C.^2 + 4*A.*B.*C*x + 4*A.*B.*C.*D*x^2;
Delta   = sqrt(temp);
y       = (-C+Delta)./(2*B.*C);
end

function y = Fn_upsilon(A,B,C,D,x,N,lambda)
temp    = C.^2 + 4*A.*B.*C*x + 4*A.*B.*C.*D*x^2;
Delta   = sqrt(temp);
f       = (-C+Delta)./(2*B.*C);
f1      = (A+2*A.*D*x)./Delta;
y1      = (A.*(1+B.*f)-A.*B.*f1*x)./((1+B.*f).^2+A.*(1+B.*f)*x)*log2(exp(1))/N;
y       = (y1-lambda)./f1;
end

function x = Fn_inv(A,B,C,D,N,lambda,upsilon,P)
x_max = P;
x_min = 0;
Delta = P;
w     = 1e-6;
Num   = 1000;
it    = 0;
while(Delta>w)
    it = it + 1;
    x  = (x_max + x_min)/2;
    temp1 = Fn(A,B,C,D,x,N);
    temp2 = Fn_lambda_upsilon(A,B,C,D,x,lambda,upsilon);
    temp1 = log(temp1);
    temp2 = log(temp2);
    if temp1 < temp2
        x_max = x;
    else
        x_min = x;
    end
    Delta = abs (temp1-temp2);
    if it == Num
        break;
    end
end
end

function y = Gn_lambda(A,B,C,D,x,N,upsilon)
y = Fn_upsilon(C,D,A,B,x,N,upsilon);
end

function y = Gn_upsilon(A,B,C,D,x,N,lambda)
y = Fn_lambda(C,D,A,B,x,N,lambda);
end

function x = gn(A,B,C,D,y)
x = fn(C,D,A,B,y);
end

function y = Gn_inv(A,B,C,D,N,lambda,upsilon,P)
y = Fn_inv(C,D,A,B,N,upsilon,lambda,P);
end










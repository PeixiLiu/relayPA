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
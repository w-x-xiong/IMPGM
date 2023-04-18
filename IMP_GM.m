function [x, fail] = IMP_GM(rangems, P0, P, gamma, anc, maxiter, eta, eps)

%Paper: A Low-Complexity Iterative Message Passing Algorithm for Robust RSS-TOA IoT Localization

%-Inputs
%anc - matrix including sensor positions 
%rangems - measured distance vector
%P - RSS vector
%P0 - initial RSS at reference distancce 1
%gamma - model parameter
%maxiter - max iteration number
%eta - threshold
%eps - a small tolerance

%-Outputs
%x - location estimate
%fail - true if fail to converge, otherwise false

[H, L] = size(anc);

x = zeros(H,1);

w = zeros(2*L,1);

sigma = zeros(H,L);
theta = zeros(L,1);

beta = zeros(L,1);

for i = 1:L
    beta(i) = 10^((P(i)-P0)/10/gamma);
end

fail = false;

%k: counter for iterations
k = 0;

while 1
    
    x_old = x;
    
    %update w, epsilon, and zeta
    for i = 1:L
        w(i) = 1/(((10*gamma/(log(10)))*(1-beta(i)*norm(x_old - anc(:,i))))^2+eps^2);
        w(i+L) = 1/((rangems(i)-norm(x_old - anc(:,i)))^2+eps^2);
        sigma(:,i) = (x_old - anc(:,i))/norm(x_old - anc(:,i));
        theta(i) = norm(x_old - anc(:,i)) - sigma(:,i)'*x_old;
    end
    
    %update x
    sum1 = zeros(H,H);
    sum2 = zeros(H,1);
    for i = 1:L
        sum1 = sum1 + (sigma(:,i)*(sigma(:,i)'))*(200*beta(i)^2*gamma^2*w(i)/log(10)^2+2*w(i+L));
        sum2 = sum2 + sigma(:,i)*(200*beta(i)^2*gamma^2*w(i)/log(10)^2+2*w(i+L))*((200*beta(i)^2*gamma^2*w(i)/log(10)^2+2*w(i+L))^(-1)*(200*beta(i)*gamma^2*w(i)/log(10)^2+2*w(i+L)*rangems(i))-theta(i));
    end
    
    if ((isnan(sum(sum(sum1)))) || (isnan(sum(sum(sum2)))))
        fail = true;
        break
    end
    
    x = (pinv(sum1))*sum2;
    
    if ((norm(x - x_old)/min(norm(x),norm(x_old))) < eta) || ((k+1) == maxiter)
        if ((k+1) == maxiter)
            fprintf('have reached the maximum iteration number\n')
            fail = true;
        end
        break
    end
    
    k = k + 1;
    
end

fprintf('It takes %d iterations to meet a termination condition\n', k)

end


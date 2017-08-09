function [P num_comb] =partition_fn(Py, k_max, N, x)
% The Forward Recursion step of the Bayesian Change Point Algorithm. P(k,j) is
% the probability of the data Y_1:j containing k change points.  num_comb is
% the total number of solutions that satisfy the constraints on the minimum 
% distance between adjacent change points - it will be used as a part of the 
% prior distribution on the number of change points in the main script.

P=zeros(k_max,N)-Inf;        % -Inf b/c starts in log form
num_comb=zeros(k_max,N);     % Counter.... used instead of n_choose_k function

%k is the number of Change Points
k=1;            % First row is different from the rest, as you add together two homogeneous segments

for j=k+1:N     % Note: Several of these terms will be -INF, due to d_min
    temp=zeros(1,j-1);
    count=0;
    for v=1:j-1
        temp(v)= Py(1,v)+Py(v+1,j);     % Changepoints occur at start of new segment, not at end of old one
        if (temp(v) > -Inf)             % If one or the other is not -Inf
            count=count+1;              % Counts the number of valid probabilities added together
        end
    end
    num_comb(k,j)=count;  
    P(k,j)=logsumlog(temp);             % Equation (3) - Marginalize over all possible placements of the change point.
end

for k=2:k_max
    for j=(k+1):N  % Note: Several of these terms will be -INF as well

    temp=zeros(1,j-1);
    temp_count=temp;            % This will take into account segments marginalized out in previous steps
    for v=1:j-1         
        temp(v) = P(k-1,v)+Py(v+1,j);
        if (temp(v) > -Inf)     %if one or the other is not -Inf, i.e. a valid solution
            temp_count(v) = num_comb(k-1,v);
        end
    end
    count=sum(temp_count);      % Total number of valid solutions
    num_comb(k,j)=count;
    P(k,j)=logsumlog(temp);     % Equation (4) - Marginalize over all possible placements of the change point.
    end
end
function P = logsumlog(P)
% Input vector P contains the log (base e) of probabilities
% Output scalar P = log(exp(P1)+exp(P2)+...+exp(Pn))
% Function ensures that numerical precision is kept in calculations involving logarithms

P = reshape(P,numel(P),1);  % Ensure a column vector
len = length(P);            % Number of elements
P = sort(P,'descend');      % Large probabilities first
while len>1;            % If at least two items remain in vector
    % Replace first two elements with sum of first two
    if (P(2)>-Inf)          % Values are sorted... so if we find -Inf then there are no more terms to add
        P=[addtwo(P(1:2));P(3:len)];
        len=len-1;          % Decrease the length
    else P=P(1);
        len=0;
    end
end % Stop when only one element left

function w = addtwo(v)
w = v(1)+log(1+exp(v(2)-v(1)));   
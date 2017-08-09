function mutualInfo = mutualInformation(counts)
% see http://stackoverflow.com/questions/6158015/mutual-information-of-matlab-matrix

  pXY = counts./sum(counts(:));
  pX = sum(pXY,2);
  pY = sum(pXY,1);

  mutualInfo = pXY.*log(pXY./(pX*pY));
  mutualInfo = sum(mutualInfo(:));

  
  
function I = MutualInformation(X,Y);
 % source: http://xaphire.de/recipes/?p=376
if (size(X,2) > 1)  % More than one predictor?
    % Axiom of information theory
    I = JointEntropy(X) + Entropy(Y) - JointEntropy([X Y]);
else
    % Axiom of information theory
    I = Entropy(X) + Entropy(Y) - JointEntropy([X Y]);
end
    
    
function H = JointEntropy(x, y, binSize)
 
    %I suggest binSize=[range(x)/std(x)*5,range(y)/std(y)*5]
    [N,C] = hist3([x y],binSize);

    hx = C{1,1};
    hy = C{1,2};

    %xyPos=meshgrid(hx,hy);

    % Normalize the area of the histogram to make it a pdf
    N = N ./ sum(sum(N));
    b=hx(2)-hx(1);
    l=hy(2)-hy(1);

    % Calculate the entropy
    indices = N ~= 0;
    H = -b*l*sum(N(indices).*log2(N(indices)));    


function H = Entropy(y,binSize)
    % Calculate the entropy for an integer value of y
 
    %I suggest: binSize=range(y)/std(y)*5;
    % Generate the histogram
    [n xout] = hist(y, binSize);
 
    % Normalize the area of the histogram to make it a pdf
    n = n / sum(n);
    b=xout(2)-xout(1);
 
    % Calculate the entropy
    indices = n ~= 0;
    H = -sum(n(indices).*log2(n(indices)).*b);

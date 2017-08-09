function Pth=OcPredFast(Plist,Pexp,hstat,Jstat,hr,Jr,h,J,J1)
%Same than OcPred, but does not estimate the empirical probability, which is the longest part. The
%pattern is given in argument. 
TempSize=size(Plist,2);

for ind=1:size(Plist,3)
        %We also store the pattern in a table (can be disabled if too slow):
    pattern(:,:)=Plist(:,:,ind);
    
    Eth=hstat'*pattern(:,1)+0.5*pattern(:,1)'*Jstat*pattern(:,1);
    for i=2:TempSize
        Eth = Eth+dot(pattern(:,i-1),J1*pattern(:,i)) + h'*pattern(:,i)+0.5*pattern(:,i)'*J*pattern(:,i) +  hr'*pattern(:,i-1)+0.5*pattern(:,i-1)'*Jr*pattern(:,i-1);
    end
    Pth(ind)=exp(Eth);
        
end


Pth=Pth*sum(Pexp)/sum(Pth);%Normalisation, faster than estimating the Z, and makes not that much difference


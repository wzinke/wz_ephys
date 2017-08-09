function [Pexp,Pth]=OcCompare(raster,Plist,hstat,Jstat,hr,Jr,h,J,J1,TempSize)
%Same type of comparison than OcPred
N=length(h);

lendata=size(raster,2);

for ind=1:size(Plist,3)
    pattern(:,:) = Plist(:,:,ind);
    test1=pattern'*raster/N;
    test2=zeros(1,lendata-TempSize+1);
    for k=1:TempSize
        test2=test2+test1(k,k:(lendata-TempSize+k));
    end
    Pexp(ind) = length(find((test2/TempSize)==1));    
	Pexp(ind) = Pexp(ind)/(lendata-TempSize+1);        
    Eth(ind)=hstat'*pattern(:,1)+0.5*pattern(:,1)'*Jstat*pattern(:,1);    
	for i=2:TempSize                
    	Eth(ind)=Eth(ind)+dot(pattern(:,i-1),J1*pattern(:,i)) + h'*pattern(:,i)+0.5*pattern(:,i)'*J*pattern(:,i) +  hr'*pattern(:,i-1)+0.5*pattern(:,i-1)'*Jr*pattern(:,i-1);
	end
	Pth(ind)=exp(Eth(ind));
end

Pth=Pth*sum(Pexp)/sum(Pth);



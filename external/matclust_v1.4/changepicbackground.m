function newpic = changepicbackground(pic,newbackground)

for i = 1:3
    temppic = pic(:,:,i);
    temppic(find(temppic ~= 0)) = newbackground(i);
    newpic(:,:,i) = temppic;
end
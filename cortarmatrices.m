function [l p]=cortarmatrices(A,sensor)
l=1;
for k=1:length(A(1,:,sensor))
    if A(1,k,sensor)==0
        l=l+1;
    else
        break
    end
end
m=1;
for k=1:length(A(1,:,sensor))-1
    if A(1,length(A(1,:,sensor))-k+1,sensor)==0 & A(1,length(A(1,:,sensor))-k,sensor)==0
        m=m+1;
    else
        break
    end
end
p=length(A(1,:,sensor))-m;
        
end
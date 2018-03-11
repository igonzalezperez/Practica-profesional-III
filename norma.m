function [An]=norma(A)
[D1 D2 D3]=size(A);
for i=1:D2
    for j=1:D3
        An(i,j)=norm(A(:,i,j));
    end
end
end

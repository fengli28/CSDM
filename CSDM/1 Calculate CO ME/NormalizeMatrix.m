function NM=NormalizeMatrix(L)
% normalization matrix L using min-max normalization

ng=size(L,1);
L1=L((L>0));
m1=max(max(L1)); % maximum value
m2=min(min(L1)); % minimum value
NM=zeros(ng,ng);
for i=1:ng
    for j=1:ng
        if L(i,j)~=0
        NM(i,j)=(L(i,j)-m2)/(m1-m2); 
        end
    end
end

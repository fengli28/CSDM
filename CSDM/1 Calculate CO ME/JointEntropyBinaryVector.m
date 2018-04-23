function JH_AB=JointEntropyBinaryVector(A,B)
% the joint entropy of two binary vectors A and B
n=size(A,2);
T=zeros(2,2); % There are 2*2=4 combinations.
% T(1,1): The number of times (0,0) appears.
% T(1,2): The number of times (0,1) appears.
% T(2,1): The number of times (1,0) appears.
% T(2,2): The number of times (1,1) appears.

for i=1:n
    if A(i)==0
        if B(i)==0
            T(1,1)=T(1,1)+1; % T(1,1): The number of times (0,0) appears.
        elseif B(i)==1
            T(1,2)=T(1,2)+1; % T(1,2): The number of times (0,1) appears.
        end
    elseif A(i)==1
        if B(i)==0
            T(2,1)=T(2,1)+1; % T(2,1): The number of times (1,0) appears.
        elseif B(i)==1
            T(2,2)=T(2,2)+1; % T(2,2): The number of times (1,1) appears.
        end
    end
end

P=T./n;
% P(1,1): The probability of (0,0).
% P(1,2): The probability of (0,1).
% P(2,1): The probability of (1,0).
% P(2,2): The probability of (1,1).

E=zeros(2,2);%
for i=1:2
    for j=1:2
        if P(i,j)~=0
            E(i,j)=-P(i,j)*log2(P(i,j));
        else
            E(i,j)=0;
        end
    end
end
JH_AB=sum(sum(E)); % the joint entropy of two binary vectors A and B


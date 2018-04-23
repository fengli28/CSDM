function H=EntropyBinaryVector(A)
% calculate the entropy of a binary vector A
% H: the entropy of A; H=-((P(1))*log(P(1))+(P(2))*log(P(2)));
% for example, A=[0 0 1 1 1]

n=length(A); % the length of vector A
N=[0 0]; % N(1): number of zeros in A; N(2): number of ones in A;
for i=1:1:n
    if A(i)==0
        N(1)=N(1)+1;
    elseif A(i)==1
        N(2)=N(2)+1;
    end
end
P=N./n; % P(1): probability of zeros in A; P(2): probability of ones in A;
E=[0 0];
for i=1:2
    if P(i)~=0
        E(i)=-P(i)*log2(P(i));
    else
        E(i)=0;
    end
end
H=sum(E); % the entropy of A


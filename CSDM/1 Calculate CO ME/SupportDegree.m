function SD=SupportDegree(A,B)
% The two vectors A and B have the same length.
% A=[0 1 1 0 1 1 1 0 0 0 0 0 1 1 1 1]
% B=[1 1 1 0 1 0 1 0 1 0 0 0 1 1 0 1]

n=size(A,2);
p=0; % the times of occurrences of A and B.
for i=1:n
    if A(i)==1 && B(i)==1
        p=p+1;
    end
end

q=0;% the times of occurrences of A or B.
for i=1:n
    if A(i)==1||B(i)==1
        q=q+1;
    end
end
SD=p/q;% the support degree of (A,B)

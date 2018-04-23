function [GenePairUC,Type]=UncertaintyCoefficentMatrix(M)
% calculate the uncertainty coefficent for each gene pair
% M: a binary muation matrix
% GenePairUC: the uncertainty coefficent and for each gene pair
% Type: the relationship between each gene pair

m=size(M,1); % number of rows (genes) 
GenePairUC=zeros(m,m);
Type=zeros(m,m);
for i=1:m
    i
    for j=i:m
        A=M(i,:); % the ith row of the matrix L
        B=M(j,:); % the jth row of the matrix L
        if sum(A)==0 || sum(B)==0 % if A or B is all zeros, let the exclusivity
            GenePairUC(i,j)=0; 
            GenePairUC(j,i)=0;
            Type(i,j)=0;
            Type(j,i)=0;
        else
            AB=or(A,B);
            S=find(AB>0); % corrected profiles: only consider the columns of A or B is ones and add one zero as a disturbance.
            s=size(S,2);
            
            A1=zeros(1,s+1); % the correct profile of A
            B1=zeros(1,s+1); % the correct profile of B
            for k=1:s
                A1(k)=A(S(k));
                B1(k)=B(S(k));
            end
            
            GenePairUC(i,j)=Compute_UC(A1,B1); % the uncertainty coefficent for (A1,B1)
            GenePairUC(j,i)=Compute_UC(B1,A1);% the uncertainty coefficent for (B1,A1)
            supp1=SupportDegree(A1,B1); % the support degree of A1 and B1
            supp2=SupportDegree(~A1,B1); % the support degree of ~A1 and B1
            supp3=SupportDegree(A1,~B1); % the support degree of A1 and ~B1

            if supp1>supp2 % the overlaps of A1 and B1 are more than that of ~A1 and B1
                Type(i,j)=1;%
            else
                Type(i,j)=2; % the overlaps of A1 and B1 are less than that of ~A1 and B1
            end
            if supp1>supp3
                Type(j,i)=1; % the overlaps of A1 and B1 are more than that of A1 and B1
            else
                Type(j,i)=2; % the overlaps of A1 and B1 are less than that of A1 and ~B1
            end
        end
    end
end
for i=1:m
    GenePairUC(i,i)=0; 
end
for i=1:m
    Type(i,i)=0; 
end

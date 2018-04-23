function MutualExclusivity_k=CalculatedMutualExclusivity(Mk)
% Mk: mutation matrix of cancer k
% MutualExclusivity_k: mutual exclusivity for each gene pair of cancer k
ng=size(Mk,1); % number of genes in cancer k
ns=size(Mk,2); % number of samples in cancer k

[GenePairUC,Type]=UncertaintyCoefficentMatrix(Mk); % calculate the uncertainty coefficent and relationship type

UC=GenePairUC.*(Type==2); % uncertainty coefficent for gene pairs 
ME=zeros(ng,ns); % mutual exclusivity for each gene pair
for i=1:ng
    i
    for j=(i+1):ng
        if UC(i,j)>0 && UC(j,i)>0 && Type(i,j)==2 && Type(j,i)==2
            ME(i,j)=(UC(i,j)+UC(j,i))/2;
            ME(j,i)=ME(i,j);
        else 
            ME(i,j)=0;
            ME(j,i)=0;
        end
    end
end

MutualExclusivity_k=NormalizeMatrix(ME);% normalization 


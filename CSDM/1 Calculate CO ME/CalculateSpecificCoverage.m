function Coverage_k=CalculateSpecificCoverage(Mk,Mall)
% Mk: mutation matrix of cancer k
% Mall: the mutation matrix of all cancer types
% Coverage_M: the specific coverage matrix of cancer k

ng=size(Mk,1); % the number of genes
nsk=size(Mk,2); % the number of samples in cancer k
GeneCM=zeros(ng,ng);% gene coverage matrix
for i=1:ng
    i
    for j=(i+1):ng
        if sum(Mk(i,:))>0 && sum(Mk(j,:))>0
        Ck=or(Mk(i,:),Mk(j,:));
        ck=sum(Ck); % the number of samples in (gi,gj) in cancer k
        Cpan=or(Mall(i,:),Mall(j,:));
        cpan=sum(Cpan); % the number of samples in (gi,gj) in pan cancer 
        GeneCM(i,j)=sqrt((ck/nsk)*(ck/cpan));
        GeneCM(j,i)=GeneCM(i,j);
        end
    end
end

Coverage_k=NormalizeMatrix(GeneCM);% normalization 



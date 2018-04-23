function [SM,SPK]=DeleteRepeatCliques(M,PK)

% M=ModulesNum; % All modules
nm=size(M,1);
nmsize=size(M,2);
E=M;
E1=M;

for i=1:nm
    i
    Mi=E(i,E(i,:)>0);
    ni=length(Mi);
    for j=(i+1):nm
        Mj=E(j,E(j,:)>0);
        nj=length(Mj);
        P=intersect(Mi,Mj);
        np=length(P);
        if ni==np % Mi is a subset of Mj
            E1(i,:)=zeros(1,nmsize); % the i th row is redundant
        elseif nj==np
            E1(j,:)=zeros(1,nmsize); % the j th row is redundant
        end
    end
end
F=find(sum(E1,2)>0);
SM=E1(F,:);
SPK=PK(F,:);



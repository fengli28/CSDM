clear all
clc

load('CO_BLCA.mat')
load('CO_BRCA.mat')
load('CO_COADREAD.mat')
load('CO_GBM.mat')
load('CO_HNSC.mat')
load('CO_KIRC.mat')
load('CO_LAML.mat')
load('CO_LUAD.mat')
load('CO_LUSC.mat')
load('CO_OV.mat')
load('CO_UCEC.mat')

load('PanGenes.mat')
ng=size(PanGenes,1);

CO_All{1,1}=CO_BLCA;
CO_All{2,1}=CO_BRCA;
CO_All{3,1}=CO_COADREAD;
CO_All{4,1}=CO_GBM;
CO_All{5,1}=CO_HNSC;
CO_All{6,1}=CO_KIRC;
CO_All{7,1}=CO_LAML;
CO_All{8,1}=CO_LUAD;
CO_All{9,1}=CO_LUSC;
CO_All{10,1}=CO_OV;
CO_All{11,1}=CO_UCEC;

nc=size(CO_All,1);

E=zeros(ng,ng);
RankFCO11={};
for i=1:nc
    RankFCO11{i,1}=E;
end

for i=1:ng
    i
    for j=(i+1):ng
        V=zeros(nc,1); % mutual exclusivity in 11 cancers and pan cancer, include pan FCO
        for k=1:nc
            V(k,1)=CO_All{k,1}(i,j); %
        end
        
        % for cancers exclude pan cancer
        % rank
        T=V;
        if sum(T)>0
            [B,I]=sort(T,'descend');
            bmax=max(B);
            Fmax=find(B==bmax); % the max value is only one cancer
            nmax=size(Fmax,1);
            
            if nmax==1
                ct=I(1,1); % which type of cancer
                RankFCO11{ct,1}(i,j)=1; % the rank for each cancer for each gene pair
                RankFCO11{ct,1}(j,i)=1;
            end
        end
    end
end
% save RankFCO11.mat RankFCO11

%% save 11 dataset
TopCO_BLCA=RankFCO11{1,1};
TopCO_BRCA=RankFCO11{2,1};
TopCO_COADREAD=RankFCO11{3,1};
TopCO_GBM=RankFCO11{4,1};
TopCO_HNSC=RankFCO11{5,1};
TopCO_KIRC=RankFCO11{6,1};
TopCO_LAML=RankFCO11{7,1};
TopCO_LUAD=RankFCO11{8,1};
TopCO_LUSC=RankFCO11{9,1};
TopCO_OV=RankFCO11{10,1};
TopCO_UCEC=RankFCO11{11,1};

save TopCO_BLCA.mat TopCO_BLCA
save TopCO_BRCA.mat TopCO_BRCA
save TopCO_COADREAD.mat TopCO_COADREAD
save TopCO_GBM.mat TopCO_GBM
save TopCO_HNSC.mat TopCO_HNSC
save TopCO_KIRC.mat TopCO_KIRC
save TopCO_LAML.mat TopCO_LAML
save TopCO_LUAD.mat TopCO_LUAD
save TopCO_LUSC.mat TopCO_LUSC
save TopCO_OV.mat TopCO_OV
save TopCO_UCEC.mat TopCO_UCEC


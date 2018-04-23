% Select the gene pairs that have the largest mutual exclusivity
% in cancer k among all cancer types
% The result TopME_BLCA is a binary matrix, and the element is 1
% if the corresponding gene pair has the largest mutual exclusivity
% among all cancer types.
clear all
clc

load('ME_BLCA.mat')
load('ME_BRCA.mat')
load('ME_COADREAD.mat')
load('ME_GBM.mat')
load('ME_HNSC.mat')
load('ME_KIRC.mat')
load('ME_LAML.mat')
load('ME_LUAD.mat')
load('ME_LUSC.mat')
load('ME_OV.mat')
load('ME_UCEC.mat')

load('PanGenes.mat')
ng=size(PanGenes,1);

ME_All{1,1}=ME_BLCA;
ME_All{2,1}=ME_BRCA;
ME_All{3,1}=ME_COADREAD;
ME_All{4,1}=ME_GBM;
ME_All{5,1}=ME_HNSC;
ME_All{6,1}=ME_KIRC;
ME_All{7,1}=ME_LAML;
ME_All{8,1}=ME_LUAD;
ME_All{9,1}=ME_LUSC;
ME_All{10,1}=ME_OV;
ME_All{11,1}=ME_UCEC;

nc=size(ME_All,1);

E=zeros(ng,ng);
RankFME11={};
for i=1:nc
    RankFME11{i,1}=E;
end

for i=1:ng
    i
    for j=(i+1):ng
        V=zeros(nc,1); % mutual exclusivity in 11 cancers and pan cancer, include pan FME
        for k=1:nc
            V(k,1)=ME_All{k,1}(i,j); %
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
                RankFME11{ct,1}(i,j)=1; % the rank for each cancer for each gene pair
                RankFME11{ct,1}(j,i)=1;
            end
        end
    end
end
% save RankFME11.mat RankFME11

%% save 11 dataset
TopME_BLCA=RankFME11{1,1};
TopME_BRCA=RankFME11{2,1};
TopME_COADREAD=RankFME11{3,1};
TopME_GBM=RankFME11{4,1};
TopME_HNSC=RankFME11{5,1};
TopME_KIRC=RankFME11{6,1};
TopME_LAML=RankFME11{7,1};
TopME_LUAD=RankFME11{8,1};
TopME_LUSC=RankFME11{9,1};
TopME_OV=RankFME11{10,1};
TopME_UCEC=RankFME11{11,1};

save TopME_BLCA.mat TopME_BLCA
save TopME_BRCA.mat TopME_BRCA
save TopME_COADREAD.mat TopME_COADREAD
save TopME_GBM.mat TopME_GBM
save TopME_HNSC.mat TopME_HNSC
save TopME_KIRC.mat TopME_KIRC
save TopME_LAML.mat TopME_LAML
save TopME_LUAD.mat TopME_LUAD
save TopME_LUSC.mat TopME_LUSC
save TopME_OV.mat TopME_OV
save TopME_UCEC.mat TopME_UCEC


clear all
clc

ImportAC=readtable('AllCancers.txt','Delimiter','tab','ReadVariableNames',true);
Genes=ImportAC{:,1};
Mall=ImportAC{:,2:end};% mutation matrix of all cancer types 

ImportM=readtable('BLCA.txt','Delimiter','tab','ReadVariableNames',true);
Mk=ImportM{:,2:end};% mutation matrix of BLCA

load('TopCO_BLCA.mat') %  rank in pan cancer of coverage
load('TopME_BLCA.mat') %  rank in pan cancer of mutual exclusivity
load CO_BLCA.mat  % the specific coverage matrix
load ME_BLCA.mat % Mutual Exclusivity matrix 
%%
SpeNetwork_Z3=SpecificNetworkConstruction(CO_BLCA,ME_BLCA,TopCO_BLCA,TopME_BLCA);

%%
[SpeModules_k,SpeMPercent_k]=SpecificModuleDetection(SpeNetwork_Z3,Mk,Mall);

%% sort driver pathways
SM=SpeModules_k;
SMIE=SpeMPercent_k; % 1: number of genes in module; 2: internal coverage; 3: external coverage
SC=zeros(size(SMIE,1),1);% specific coverage
for i=1:size(SMIE,1)
    SC(i,1)=sqrt(SMIE(i,2)*SMIE(i,3)); % specific coverage
end
[F,IF]=sort(SC(:,1),'descend');% rank by specific coverage
Modules_k=SM(IF,:); % specific modules
SC_k=SC(IF,:); % modules' number of genes and specific coverage
InEx=SMIE(IF,:);
% save Modules_k.mat Modules_k
% save PercentM_k.mat PercentM_k

txtoutModules(Genes,Modules_k,SC_k,InEx)

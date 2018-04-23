function [SpeModules_k,SpeMPercent_k]=SpecificModuleDetection(SpeNetwork,Mk,Mall)
% SpeNetwork: a binary matrix to represent the specific matrix
% Mk: the mutation matrix of cancer k
% Mall: the mutation matrix of all cancer types
% SpeModules_k: the genes involved in the modules
% SpeMPercent_k: three columns, 1: number of genes in the module;
% 2: the external coverage of the module; 3: the internal coverage of the
% module.

% Mk=MutBLCA; % mutation matrix of cancer k
% Mall=PanMut; % mutation matrix in pan cancer

ng=size(Mk,1); % the number of genes
nsk=size(Mk,2); % the number of samples in cancer k

SN=SpeNetwork;
for i=1:ng
    SN(i,i)=1; % 
end

%% specific module detection
MEModulesAll={}; % the set of specific modules
PerKpan=zeros(ng,3);
for i=1:ng
    i
    
    Module_i=i; % the module i including gene i as a seed
    ProfileK_i=Mk(i,:); % the profile of  gene i in cancer k
    ProfilePan_i=Mall(i,:); % the profile of  gene i in pan cancer
    cover_ki=sum(ProfileK_i); % number of samples with mutated gene i in cancer k
    cover_pani=sum(ProfilePan_i); % number of samples with mutated gene i in all cancers
    percent_kpan=cover_ki/cover_pani; % the percentage of samples in cancer i
    
    NeighborsK=find(SN(i,:)>0); % the neighbors of  gene i in cancer k
    NeighborsK=setdiff(NeighborsK,Module_i); % the neighbors of gene i exclude i
    
    if isempty(NeighborsK) % gene i is an isolated point
        MEModulesAll{i,1}=Module_i;
        continue
    end
    
    A=ProfileK_i;
    % k=module_size; 
    k=length(NeighborsK); % the maximum size of modules
    while k>0
        nb=size(NeighborsK,2); %the number of neighbors of gene i
        ProfileK_Nb=Mk(NeighborsK,:); % the profiles of all neighbors of gene i in cancer k
        ProfilePan_Nb=Mall(NeighborsK,:); % the profiles of all neighbors of gene i in pan cancer
        
        PercentNb=zeros(nb,1);% percentage of samples with each neighbor in pan cancer
        for j=1:nb
            Nj=NeighborsK(j); % the j th neighbor of module i
            Pk_Nj=ProfileK_Nb(j,:); % the  mutated profile of j th neighbor
            Ppan_Nj=ProfilePan_Nb(j,:); % the  mutated profile of j th neighbor
            Ck=Mk(Module_i,:);
            Cpan=Mall(Module_i,:);
            SNjModiK=[Pk_Nj;Ck];
            SNjModiPan=[Ppan_Nj;Cpan];
            Sk=length(find(sum(SNjModiK,1)>0));
            Span=length(find(sum(SNjModiPan,1)>0));
            
            p=Sk/Span;
            PercentNb(j,1)=p;
        end
        
        max_per=max(PercentNb); % decide the percentage of samples if increase
        if percent_kpan>max_per
            break
        else
            Select_j=NeighborsK(find(PercentNb==max_per)); % select the neighbor with highest  score
            
            Module_i=union(Module_i,Select_j(1)); % add the j th neighbor into module_i
            Profile_j=Mk(Select_j(1),:); % the profile of selected j th domain
            A=or(Profile_j,A);
            NewNeighbors=find(sum(SN(Module_i,:),1)==size(Module_i,2));
            NeighborsK=setdiff(NewNeighbors,Module_i); % the common neighbors for all genes in module i
            
            if isempty(NeighborsK)
                break
            end
        end
        k=k-1;
    end
    i
    Ck=Mk(Module_i,:);
    Cpan=Mall(Module_i,:);
    Sk=length(find(sum(Ck,1)>0));
    Span=length(find(sum(Cpan,1)>0));
    PerKpan(i,1)=length(Module_i);%  number of genes in module
    PerKpan(i,2)=Sk/nsk; % internal coverage
    PerKpan(i,3)=Sk/Span; % external coverage
    MEModulesAll{i,1}=Module_i; % module involved gene i
end

F3=find(PerKpan(:,1)>2); % filter the modules less than 2 genes
MEModules3=MEModulesAll(F3,:);
PerKpan3=PerKpan(F3,:);
[B,I]=sort(PerKpan3(:,3),'descend');
PerKpan3Sort=PerKpan3(I,:);
SpeModules=MEModules3(I,:);

%% merge modules
nm=size(SpeModules,1);
nmsize=max(PerKpan3Sort(:,1));
ModulesNum=zeros(nm,nmsize);
for i=1:nm
    M1=SpeModules{i,1};
    ModulesNum(i,1:length(M1))=M1;
end
[SM,SPK]=DeleteRepeatCliques(ModulesNum,PerKpan3Sort); % delete recurrent modules

SpeModules_k=SM; % the specific driver modules for cancer k
SpeMPercent_k=SPK; % 





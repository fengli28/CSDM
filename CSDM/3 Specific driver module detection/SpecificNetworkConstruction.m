function SpeNetwork_Z3=SpecificNetworkConstruction(CO,ME,TopCO,TopME)
% CO: the specific coverage matrix for cancer k
% ME: the mutual exclusivity matrix for cancer k
% RankCO: the gene pairs of maximum specific coverage across all cancer 
% types in % cancer k
% RankME: the gene pairs of maximum mutual exclusivity across all cancer 
% types in cancer k
% SpeNetwork_Z3: the unweighted specific network

ng=size(CO,1);
MaxR=(TopCO==1).*(TopME==1); % the gene pairs of both maximum specific coverage and mutual exclusivity  
COMax=CO.*MaxR; % the coverage matrix
MEMax=ME.*MaxR; % Mutual Exclusivity matrix 

%% weights for gene pairs, harmonic mean of specific coverage and mutual exclusivity
HarmonicCOME=zeros(ng,ng); 
for i=1:ng
    i
    for j=1:ng
        if MEMax(i,j)~=0 && COMax(i,j)~=0 
            HarmonicCOME(i,j)=2/(1/MEMax(i,j)+1/COMax(i,j));
        else 
            HarmonicCOME(i,j)=0;
        end
    end
end

HarmonicCMNom=HarmonicCOME/max(max(HarmonicCOME)); % normalization
clear CM
clear MEM
clear HarmonicCOME

%% cancer specific network
TrilWeight=tril(HarmonicCMNom);
V=TrilWeight(TrilWeight>0);
[Z,mu,sigma]=zscore(V);
ZscoreNet=zeros(ng,ng);
for i=1:ng
    for j=1:ng
        if HarmonicCMNom(i,j)>0
            ZscoreNet(i,j)=(HarmonicCMNom(i,j)-mu)/sigma;
        end
    end
end

%% z-score larger than 3
ZscoreNet3=HarmonicCMNom.*(ZscoreNet>3);
SpeNetwork_Z3=ZscoreNet3>0; % the relationship between genes




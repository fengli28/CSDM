clear all
clc
% calculate the specific coverage matrix and mutual exclusivity matrix
% BLCA as an example

ImportAC=readtable('AllCancers.txt','Delimiter','tab','ReadVariableNames',true);
Genes=ImportAC{:,1};
Mall=ImportAC{:,2:end};% mutation matrix of all cancer types 

ImportM=readtable('BLCA.txt','Delimiter','tab','ReadVariableNames',true);
Mk=ImportM{:,2:end};% mutation matrix of BLCA

%% calculate the specific coverage matrix
CO_BLCA=CalculateSpecificCoverage(Mk,Mall);
save CO_BLCA.mat CO_BLCA

%% calculate the mutual exclusivity matrix
ME_BLCA=CalculatedMutualExclusivity(Mk);
save ME_BLCA.mat ME_BLCA




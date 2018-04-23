function txtoutModules(Genes,Modules, SC,InEx)
% Genes: gene symbols
% Modules: modules with gene ids
% PercentM: 

fid=fopen('Modules of cancer k.txt','w+');
ModuleNames={};
fprintf(fid,'Number of genes in module\tSpecific coverage\tInternal coverage\tExternal coverage\tGenes in module\n');
for i=1:size(Modules,1)
    fprintf(fid,'%d\t%.3f\t%.3f\t%.3f\t', InEx(i,1),SC(i,1),InEx(i,2),InEx(i,3));
    for j=1:size(Modules,2)
        g=Modules(i,j);
        if g>0
            ModuleNames{i,j}=Genes{g,1};
            fprintf(fid,'%s\t', ModuleNames{i,j});
        end
    end
     fprintf(fid,'\n');
end
fclose(fid);

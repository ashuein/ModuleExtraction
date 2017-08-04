% Example

% Dynamic_network_generation is from TSOCD source code, refer to "Detecting temporal protein complexes from dynamic protein-protein interaction networks"
% Data_set = dynamic_network_generation('DIP_static_ppi.txt', 'DIP_expression_data.txt', delta);
addpath('../ModuleExtraction')
load('data/Data_set.mat');
GENE=Data_set.Protein;
T=12;
W=Data_set.DPIN{1};
z=Data_set.Expression(:,1);
for t = 2 : T
    W = W+Data_set.DPIN{t};
    z = z+Data_set.Expression(:,t);
end
W(find(W<(0.5*T)))=0;
[s, subset] = ModuleExtract(W,z',1,3,10,GENE,'PPI_complexes.txt');

for t = 1 : T
    W = Data_set.DPIN{t};
    z = Data_set.Expression(:,t);
    [s, subset] = ModuleExtract(W,z',1,3,10,GENE,['PPI_complexes' num2str(t) '.txt']);
end
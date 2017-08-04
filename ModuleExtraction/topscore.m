%% FUNCTION TOPSCORE
%   calculated the highest-scoring connected component with given nodes set
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   randomscore: n*2 matrix, the i-th row are the mean and sd of random
%   nodeset: given nodes set

%% OUTPUT
%   module score and component nodes set
%
%% LICENSE
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2015 - 2016 Dong Li
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via donggeat@gmail.com
%
%   Last modified on Aug 04, 2017.
%
%% Related papers
%  [1] Discovering regulatory and signalling circuits in molecular interaction networks. Trey Ideker et al, Bioinformatics 2002
%  [2] Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. Dong Li et al, BMC Genomics 2017
%  [3] Extracting active modules from multilayer PPI network: a continuous optimization approach. Dong Li et al, in prearing
function [s,subset] = topscore(G,nodeset,minmodulesize,maxmodulesize,GENE,filename)
global StackC;
global ptr;

if nargin < 2
    error('\n Inputs: G, nodeset should be specified!\n');
end

s = -inf;
[Cnew,Lnew] = conncomp(G(nodeset,nodeset));
labels = unique(Lnew);

for i = 1:length(labels)
    nodeList = nodeset(find(Lnew==labels(i)));
    k=length(nodeList);
    
    if k > s
        s = k;
        subset = nodeList;
    end
    
    if k > maxmodulesize
        if ptr==0
            ptr = 1;
        end
        StackC{ptr,1}= nodeList;
        ptr = ptr + 1;
    elseif k >=minmodulesize
        fid = fopen(filename, 'a+');
        for j = 1:length(nodeList)
            fprintf(fid, '%s\t', GENE{nodeList(j)});
        end
        fprintf(fid, '\n');
        fclose(fid);
    end
end

end

%http://stackoverflow.com/questions/16883367/how-to-find-connected-components-in-matlab
function [S,C] = conncomp(G)
  [p,q,r] = dmperm(G'+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end

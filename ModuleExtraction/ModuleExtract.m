%% FUNCTION ModuleExtract
%   Continuous optimization for active modules identification in molecular interaction networks
%
%% INPUT
%   G: n*n adjacency matrix of the network
%   array_basic_z: n*1 vecor, z-score of the nodes
%   minmodulesize: minima module size
%   maxmodulesize: maximal module size
%   GENE: n*1 cell, gene or protein names
%   filename: for saving result, each line as a module/complex
%% OUTPUT
%   Module score and nodes for the last module
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
%   Copyright (C) 2015 - 2017 Dong Li
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via donggeat@gmail.com
%
%   Last modified on Aug 04, 2017.
%
%% Related papers
%  Extracting active modules from multilayer PPI network: a continuous optimization approach. Dong Li et al, in prearing


function [s, subset] = ModuleExtract(G,array_basic_z,lambda,minmodulesize,maxmodulesize,GENE,filename)
maxiter=1000;
n = size(G,1);
subset=1:n;

% Option 2, Binary extraction
% If we want to take top-10 modules, start to tract  
% when the rest size is 10*100

N = 100;
condition = true;
global StackC;
global ptr;

StackC = cell(N,1);
ptr = 1;

while condition
        if (length(subset) >= maxmodulesize)
            Wp = G(subset,subset);
            zp = array_basic_z(subset)';
            [x, func] = MEOpt(Wp,zp,lambda,maxiter);
            tmpset = find(x<0);
            % Push it into a stack
            
            StackC{ptr,1}= subset(tmpset);
            ptr = ptr + 1;
            
            % Further partition
            tmpset = find(x>0);
            [s,subset] = topscore(G,subset(tmpset),minmodulesize,maxmodulesize,GENE,filename);
            %subset = subset(tmpset);
        elseif( length(subset) <= maxmodulesize & length(subset) >= minmodulesize)
            % save it
            fid = fopen(filename, 'a+');
            for j = 1:length(subset)
                fprintf(fid, '%s\t', GENE{subset(j)});
            end
            fprintf(fid, '\n');
            fclose(fid);
            condition = false;
        else
            disp(['discard ' num2str(length(subset))]);
%             G(subset,:)=[];
%             G(:,subset)=[];
%             array_basic_z(subset)=[];
            condition = false;
        end
end


ptr = ptr - 1;

while ptr > 0
subset=StackC{ptr,1};
StackC{ptr,1}=[];
ptr = ptr - 1;
stackcondition = true;
 while stackcondition
        if( length(subset) >= maxmodulesize)
            % Push it into a stack
            Wp = G(subset,subset);
            zp = array_basic_z(subset)';
            [x, func] = MEOpt(Wp,zp,lambda,maxiter);
            tmpset = find(x<0);
            % Push it into a stack
            if(length(tmpset) >=minmodulesize)
                ptr = ptr+1;
                StackC{ptr,1}= subset(tmpset);
            end
            
            
            % Further partition
            tmpset = find(x>0);
            [s,subset] = topscore(G,subset(tmpset),minmodulesize,maxmodulesize,GENE,filename);
            %subset = subset(tmpset);
        elseif( length(subset) <= maxmodulesize & length(subset) >= minmodulesize)
            % save it
            fid = fopen(filename, 'a+');
            for j = 1:length(subset)
                fprintf(fid, '%s\t', GENE{subset(j)});
            end
            fprintf(fid, '\n');
            fclose(fid);
            stackcondition = false;
%         elseif( length(subset) <= minmodulesize & length(subset) > 0)
%             disp(['discard ' num2str(length(subset))]);
%             G(subset,:)=[];
%             G(:,subset)=[];
%             array_basic_z(subset)=[];
%             %ptr = ptr-1;
%             stackcondition = false;
        else
            %ptr = ptr-1;
            stackcondition = false;
        end
end
end
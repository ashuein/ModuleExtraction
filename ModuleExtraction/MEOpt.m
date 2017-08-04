%% Function
%   Optimization for module extraction
%
%% INPUT
%   W: n*n adjacency matrix of the network
%   z: n*1 vecor, z-score of the nodes
%   lambda: balance parameter of two parts
%   maxiter: maximal iteration for searching proper a
%% OUTPUT
%   Module membership and objective values
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
%   Copyright (C) 2017 Dong Li
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Dong Li via donggeat@gmail.com
%
%   Last modified on Aug 04, 2017.
%
%% Related papers
%  Extracting active modules from multilayer PPI network: a continuous optimization approach. Dong Li et al, in prearing

function [x, func] = MEOpt(W,z,lambda,maxiter)

    n = size(W,1);
    
        
    k = sum(W,1);
    D = diag(k);
    L = D-W;
    
    
    epsilon = 1e-9;
    
    
    eta1 = ones(n,1);
    
    x0 = randn(n,1);
    x0 = sum(eta1)*x0-sum(x0)*eta1;
    x0 = x0/sqrt(sum(x0.^2))*sqrt(n);
    
    x = x0;
    f_x = 0.5*(x'*L*x)-lambda*(z'*x);
    grad = L*x-lambda*z;
    func = zeros(maxiter,1);
    alpha = 0.01; 
    for iteration =1:maxiter
        func(iteration) = f_x;

        % search step size alpha, not working well
         sigma = 0.01;
         beta = 0.1;
        % 
        for inner_iter=1:20,
            xp = x-alpha*grad;
            xp = sum(eta1)*xp-sum(xp)*eta1;
            xp = xp/sqrt(sum(xp.^2))*sqrt(n);
        
            %0.5*(xp'*L*xp)/n-lambda*(z'*xp)- 0.5*(x'*L*x)/n+lambda*(z'*x)-sigma*grad'*dx;          
            dx = xp-x;
            suff_decr = 0.5*(dx'*L*dx)-lambda*(z'*dx)-sigma*grad'*dx < 0;

            if inner_iter==1,
                decr_alpha = ~suff_decr;
            end
            if decr_alpha,
                if suff_decr,
                    break;
                else
                    alpha = alpha * beta;
                end
            else
                if ~suff_decr,
                    break;
                else
                    alpha = alpha/beta;
                end
            end
                
        end
        
        x_cand = x-alpha*grad;    
        % orthogonalization: sum(\|x\|)=0, construct orthogonal vector
        x_cand = sum(eta1)*x_cand-sum(x_cand)*eta1;
        % normalization: \|x\|_2^2=n
        x_cand = x_cand/sqrt(sum(x_cand.^2))*sqrt(n);

        %disp(sum(abs(x_cand-x).^2)^(1/2));
        if(sum(abs(x_cand-x).^2)^(1/2) < epsilon)
            break;
        end
        x=x_cand;
        grad = L*x-lambda*z;
        f_x = 0.5*(x'*L*x)-lambda*(z'*x);
    end
    func=func(1:iteration);
%end
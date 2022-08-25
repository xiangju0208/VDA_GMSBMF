function [M_recovery] = A_VDA_GMSBMF(matDV, Wdd, Wvv, gm,w, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter)
% % % k= floor( min(size(matDV)) * 1)
    w1=w; w2=w;
    [Gdd,Gvv] = getGIPSim_IN(matDV, gm, gm,  0, 0   );
    Sdd=w1*Gdd+(1-w1)*Wdd;
    Svv=w2*Gvv+(1-w2)*Wvv;                       
%                     lambda1 = 0.9; lambda2 = 0.8; lambda3 = 0.8;  tol1 = 2*1e-30;tol2 = 1*1e-40;      
% % %                     lambda1 = 0.1; lambda2 = 0.1; lambda3 = lambda2; k = floor(min(Nd, Nv) * 0.7); tol1 = 2*1e-3;tol2 = 1*1e-4;      
    [U, V, iter] = A_MSBMF_IN(matDV, Sdd, Svv, lambda1, lambda2, lambda3, k , tol1, tol2, maxiter); 
    M_recovery = U * V'; 

end

function [X, Y, iter] = A_MSBMF_IN(M, D, R, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter)
% MSBMF
% Usage: [X, Y, iter] = MSBMF(M, D, R, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter)
%
% Inputs:
%        M                  - the target matrix with only known entries and the unobserved entries are 0.
%        D                  - disease similarity matrix
%        R                  - drug similarity matrix
%        lambda1,2,3        - parameters needed to give.
%        k                  - the latent dimension of matrix factorization.
%        tol1, tol2         - tolerance of termination conditions.
%        maxiter            - maximum number of iterations.
%
% Outputs:
%        X, Y               - two latent low-rank matrices of the completed matrix.
%        iter               - the number of iterations.
%
% Written by: Mengyun Yang
% Email: mengyunyang@csu.edu.cn
% Created: December 16, 2019

rand('state', 2019); % fix random seed
omega = double(M ~= 0);
omega_ = ones(size(omega)) - omega;
U = rand(size(M, 1), k);
V = rand(size(M, 2), k);
P = rand(size(D, 2), k);
Q = rand(size(R, 2), k);
X = U;
Y = V;
Z = M;
W1 = zeros(size(U));
W2 = zeros(size(V));
XY = M;

rho = 1.05;
mu = 1e-4;
max_mu = 1e20;

stop1 = 1;
stop2 = 1;

for i = 1: maxiter
    U = (Z * V + lambda2 * D * P - W1 + mu * X) * inv(V' * V + lambda2 * P' * P + (lambda1 + mu) * eye(k));
    
    V = (Z' * U + lambda2 * R * Q - W2 + mu * Y) * inv(U' * U + lambda2 * Q' * Q + (lambda1 + mu) * eye(k));
    
    P = (lambda2 * D' * U) * inv(lambda2 * U' * U + lambda3 * eye(k));
    
    Q = (lambda2 * R' * V) * inv(lambda2 * V' * V + lambda3 * eye(k));
    
    X = U + (1 / mu) * W1;
    X(X < 0) = 0;
    
    Y = V + (1 / mu) * W2;
    Y(Y < 0) = 0;
    
    Z = M .* omega + (U * V') .* omega_;
    
    W1 = W1 + mu * (U - X);
    
    W2 = W2 + mu * (V - Y);
    
    stop1_0 = stop1;
    stop1 = norm(X * Y' - XY, 'fro') / norm(XY, 'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    XY = X * Y';
    if stop1 < tol1 && stop2 < tol2
        iter = i;
        break
    else
        iter = i;
        mu = min(mu * rho, max_mu);
    end
    
end

end


 %%
function [ kd,km ] = getGIPSim_IN(Adm_interaction, gamma0_d, gamma0_m,  AvoidIsolatedNodes, RemoveNonoverlapPairs   )
% Xiang  2019-11-16 
% Ref: van Laarhoven, Twan, Sander B. Nabuurs, and Elena Marchiori.
% "Gaussian interaction profile kernels for predicting drug?¡ìCtarget interaction." 
% Bioinformatics 27, no. 21 (2011): 3036-3043.
%interaction: relation matrix between disease and miRNA, row:disease    column:miRNA  
% if ~exist('gamma0_d','var') || isempty(gamma0_d)
%    gamma0_d = 1 ;  
% end
% if ~exist('gamma0_m','var') || isempty(gamma0_m)
%    gamma0_m = 1 ;  
% end
% if ~exist('AvoidIsolatedNodes','var') || isempty(AvoidIsolatedNodes)
%    AvoidIsolatedNodes = 1 ;  
% end
% % % % % % % % % % % % % % % % % % % 
if isempty( gamma0_d  ) && isempty( gamma0_m  )
    error( 'both gamma0_d and gamma0_m are empty. No output.'   );
end
[nd_all, nm_all] = size( Adm_interaction  ); 
if AvoidIsolatedNodes
    nodes_d = sum(Adm_interaction,2)~=0 ; 
    nodes_m = sum(Adm_interaction,1)~=0 ; 
    Adm=double( Adm_interaction( nodes_d , nodes_m ) );
else
    Adm=double( Adm_interaction );
end
%
if ~exist( 'RemoveNonoverlapPairs','var'  ) || isempty( RemoveNonoverlapPairs  )
    RemoveNonoverlapPairs = true ;
end

% % % % % % % % % % % % % % % % % % % 
[nd, nm] = size( Adm  );  
SumOfSquares = sum(  Adm(:).^2 ); 
% % calculate gamad for Gaussian kernel calculation
if ~isempty( gamma0_d  )
    gamma_d = gamma0_d./( ( SumOfSquares ) / nd );
    %calculate Gaussian kernel for the similarity between disease: kd
    D   = Adm*Adm';
    dd  = diag( D  ); 
    kd2 = exp(-gamma_d*( dd+(dd'-2*D)   )  );  
    if RemoveNonoverlapPairs
        kd2(~D) = 0 ;        
    end 
    if AvoidIsolatedNodes  
        kd = zeros( nd_all  );  
        kd(  nodes_d , nodes_d ) =  kd2 ;       
    else
        kd = kd2;     
    end
    %
    D = [];
    kd2 = []; 
end

% % calculate gamam for Gaussian kernel calculation
if ~isempty( gamma0_m  )
    gamma_m = gamma0_m./( ( SumOfSquares ) / nm );
    %calculate Gaussian kernel for the similarity between miRNA: km
    E=Adm'*Adm;
    mm  = diag( E  ); 
    km2 = exp(-gamma_m*( mm+(mm'-2*E)   )  );  
    %
    if RemoveNonoverlapPairs
        km2(~E) = 0 ;        
    end    
    % 
    if AvoidIsolatedNodes   
        km = zeros( nm_all  );  
        km(  nodes_m , nodes_m ) =  km2 ;      
    else 
        km = km2;     
    end
    E = [];
    km2 = [];     
end

end

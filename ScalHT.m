function [X_rec,err,time] = ...
    ScalHT(Xs,X_star,Omega,p,n,r,s,maxit,tol,mode,eta)
% ScaledGD for Hankel tensor completion (ScalHT)
%Inputs  
%    Xs:measurements.
%X_star: true matrix.
% Omega: sampling matrix  
%     p: sampling ratio
%     r: model order of the spectrally sparse signal.
%     s: the dimension of channels/multiple measurements.
% maxit: maximum number of iterations.
%   tol: exit tolerance
%  mode; run mode. 0,test run time mode. exit when mamximum iteration is reached or diverge. 
%        1, phase trans test mode: addtional exit conditions as: exit when converge/diverge/recovery error is reached.  
%
%Outputs
% X_rec: recoverd matrix
%   err:  relative error.
%   time: run time versur iteration

% clear; 
% close all; clc;
% profile on 
% profile clear
% n = 1023;
n1 = round((n+1)/2);
n2 = n+1-n1;

T = maxit;

init_mode = 0; % 0 deontes spectral init; 1 deontes random init
[s_ind ,n_ind]= find(Omega==1);
thresh_up = 1e3; 

timer = zeros(T,1);
errors = zeros(T,1);
V_store = cell(T,1);
B_store = cell(T,1);

w = zeros(n, 1);
succ_err = 1e-8;

for k = 1:n
    w(k) = min([k, n1+n2-k, n1, n2]); %length of skew-diagonals
end

   Ys = Xs.*sqrt(w')/p;

    %% Sequential Spectral initialization    
    U0 = cell(2, 1);
    if init_mode ==0
    [L0, ~, ~] = svds(@(v,tflag)blockHankelVecMul(Xs/p, v, tflag, n, n1, n2, s,1), [n1, s*n2], r); % mode 1 block Hankel from mode 1 matricization 
%     S0 = ttm(Xs_tensor/p, U0, 't');% this is conjugate transpose.
    [R0, ~, ~] = svds(@(v,tflag)blockHankelVecMul(Xs/p, v, tflag, n, n1, n2, s,2), [n2, s*n1], r); % mode 2 block Hankel from mode 2 matricization

    interm = reshape(HankelTmul(Ys,L0,2,n1,n2,n,w),[n2 r s]); % HankelT multiplication, and reshape the result n2*rs matrix into n2*r*s tensor
    interm = permute (interm,[2 1 3]); % permute this tensor's order, into a r*n2*s tensor
    interm_3 = double(tenmat(tensor(interm),3)); % mode-3 matricization
    [V0, ~,~]=svds(interm_3,r);
    % fast computation for S0 = (L',R',V') dot Ys
    W0 = conj(conv_lr(L0,R0,n,w)); % construction of W0 via r^2 fast convolution. 
    S0_3 = (V0'* Ys)* W0 ; % mode_3 matricization result of S0
    S0 = tensor(reshape(S0_3.',[r r r]));   
    else   % random init
    eta = 0.1;      
    L0 = randn(n1, r)/sqrt(n1);
    R0 = randn(n2, r)/sqrt(n2);
    V0 = randn(s, r)/sqrt(s);   
    S0_seed = randn(r, r, r);
    S0_seed = tensor(S0_seed/norm(S0_seed(:)));
%     S0 = 0.5*S0_seed*norm(Xs.*sqrt(w'),'fro')/sqrt(p);
    S0 = 1*S0_seed*norm(Xs.*sqrt(w'),'fro')/sqrt(p);
    end
    %% ScaledGD-Hankel tensor
%     U = U0;
    L = L0;
    R = R0;
    V = V0;
    S = S0;
    Ys = sparse (Ys); % Here, Ys = p^{-1}P_Omega(Y) 

tic
        % prepare the efficient representations for the 1st iteration
        W_3 = conj(conv_lr(L,R,n,w)); % construction of W via r^2 fast convolution. 
        B = W_3*double(tenmat(S,3))'; % O(nr^3), n*r  matrix  
    for t = 1:T
        Z_Omega = sum(V(s_ind,:).*conj(B(n_ind,:)),2)/p; %  O(mr),P_Omega(V*B^H)/p, only m entries of V*B' are computed£¬
        % it is with lower computation theoretically, while it runs longer than the matrix multiplication of the same complexity in matlab.
        % maybe improved or optimized specifically in matlab.  serverlal vector index and hadmard prodoct complxity (constants)         
        M_sparse = sparse(s_ind,n_ind,Z_Omega,s,n)-Ys;  
        L_sq = L'*L;
        R_sq = R'*R; % O(nr^2)
        V_sq = V'*V; % O(sr^2)
           
        M_tilde1 = V'* M_sparse; % r*n matrix, O(mr)
        M_tilde2 = (V'*V)*B'; % r*n matrix, O((s+n)r^2)
        E_hat = M_tilde1 -M_tilde2;  %r*n matrix
        
        L_update = HankelTmul(E_hat,R,1,n1,n2,n,w)*double(tenmat(S,1))'; % HankelT multiplication mode=1,O(nr^2log n+nr^3)
        L_precond = double(tenmat(ttm(S, {R_sq,V_sq}, [2 3]),1))*double(tenmat(S,1))'; % L_precond=L_breve'*L_breve, O(r^4)
        L_plus = L - eta * (L_update)/L_precond - eta*L;
        
        R_update = HankelTmul(E_hat,L,2,n1,n2,n,w)*double(tenmat(S,2))'; % HankelT multiplication mode=2,O(nr^2log n+nr^3)
        R_precond = double(tenmat(ttm(S, {L_sq,V_sq}, [1 3]),2))*double(tenmat(S,2))'; % R_precond=R_breve'*R_breve, O(r^4)
        R_plus = R - eta * (R_update)/R_precond - eta * R;       
        
        V_update1 = M_sparse*B; % O(mr)          
        V_update2 = V*(B'*B); %  O((s+n)r^2)      
        V_precond = double(tenmat(ttm(S, {L_sq,R_sq}, [1 2]),3))*double(tenmat(S,3))'; % V_precond=V_breve'*V_breve, O(r^4)       
        V_plus = V - eta*(V_update1-V_update2)/V_precond-eta*V;
        
        S_update_3 = E_hat*W_3; % r*r^2 matrices, M_3(S_update),3rd dim, O(nr^3)
        S_update = tensor(reshape (S_update_3.',[r r r]));
        S_plus = S -eta * ttm(S_update,{(L_sq)^(-1) ,(R_sq)^(-1),(V_sq)^(-1)}) -eta * S; % O(r^4)
        
        % prepare the efficient representations for the next iteration/output
        W_3 = conj(conv_lr(L_plus,R_plus,n,w)); % construction of W via r^2 fast convolution. 
        B = W_3*double(tenmat(S_plus,3))'; % O(nr^3), n*r  matrix       
        timer(t) = toc;    
        time = timer(1:t);  
              
        if (mode ==1) % phase transition test mode
        ratio1 = norm(L_plus-L,'fro')/norm(L,'fro');
        ratio2 = norm(R_plus-R,'fro')/norm(R,'fro');
        ratio3 = norm(V_plus-V,'fro')/norm(V,'fro');
        ratio4 = norm(S_plus-S)/norm(S);
        ratio = max([ratio1,ratio2,ratio3,ratio4]);      
%         fprintf('iter: %d, ratio %8f:. \n',t,ratio);
        Z = V*B'; % O(nsr)       
        X = Z.*(1./sqrt(w'));
        X_rec = X;
        error = norm(X - X_star,'fro')/norm(X_star,'fro');
        errors(t) = error;
        err = errors(1:t);
            if ~isfinite(error) || error > thresh_up || error < succ_err || ratio <= tol
                break;
            end
        end
        L = L_plus;
        R = R_plus;
        V = V_plus;
        S = S_plus;
        if (mode ==0) % test run time mode
        V_store{t} = V;
        B_store{t} = B;
        end
        if norm(S) > 1e6
            fprintf('ScalHT Diverge. \n')
        break;
        end
    end
    %% we calculate the relative error beyond the trjectories
     if mode == 0
         for j = 1:t
            Z = V_store{j}*B_store{j}'; % O(nsr)       
            X = Z.*(1./sqrt(w'));
            error = norm(X - X_star,'fro')/norm(X_star,'fro');
            errors(j) = error;
            err = errors(1:j);
         end
         X_rec = X;
%          err = error;
    end
    
end

function  [W3]=conv_lr(L,R,n,w)
    % return W3: Mode-3 matricization of W, an n*r^2 matrix. 
    % convolution of colums of L (n1*r matrix) and R (n2*r matrix)
    % via a kronector manner in the row dim, which are r^2 convolution in
    % toal
    % return a n*r^2 matrix W3
    nn = 2^nextpow2(n); % padded to a length of 2^k .
    L_fft = fft(L, nn);
    R_fft = fft(R, nn);
%     L_col = size(L_fft,2);
%     L_rep = repmat(L_fft,1,L_col); %nn * r^2 matrix
%     R_fft = fft(R, nn);
%     R_col = size(R_fft,2);
%     R_kron = kron(R_fft,ones(1,R_col)); % nn * r^2 matrix
%     W3 = ifft(L_rep.* R_kron);
    W = khatrirao(R_fft.',L_fft.'); % Khatri-Rao product
    W3 = ifft(W.');
    W3 = W3(1:n,:).*(1./sqrt(w));
end

function [result]=HankelTmul(M,F,mode,n1,n2,n,w)
% input: M, r * n matrix / s * n matrix
%mode 1: return  M_1 (G(M)*_2 R'), n1 * r^2 / n1 * (rs) matrix, 2sr dim's order is
%vec([r s])
%mode 2 :return  M_2 (G(M)*_1 L')  n2 * r^2 / n2 * (rs) matrix
    M = M.';
    M = M.*(1./sqrt(w));
    nn = 2^nextpow2(n); % padded to a length of 2^k .
    M_fft = fft(M,nn);
    F_fft =  fft(flip(conj(F)),nn);
    inter_result = khatrirao(M_fft.',F_fft.'); % Khatri-Rao product, left kronecker right
    result = ifft(inter_result.');   
%     M_col = size (M,2);
%     M_fft = fft(M,nn);
%     M_kron = kron(M_fft,ones(1,M_col));    
%     F_col =size (F,2);
%     F_fft =  fft(flip(conj(F)),nn);
%     F_rep = repmat(F_fft,1,F_col);
%     result = ifft(F_rep.*M_kron);
    if mode==1
    result = result(n2:n,:);    
    else 
    result = result(n1:n,:);    
    end
end

%     
function z = blockHankelVecMul(M, v, tflag, nd, n1, n2,s,mode)
% block Hankel multiplication Z=M_1(H[M])v if tflag='notransp'; z=M_1(H[M])'v if tflag='transp'
% mode 1, matriczation 1 , is a n1*n2 s block Hankel matrix; mode 2, is a n2*n1s block Hankel matrix. 
    if mode ==1
        na = n1;
        nb = n2;
    else 
       na = n2;
       nb = n1;
    end
    n = 2^nextpow2(nd);
    z = zeros(na,1);
    if strcmp(tflag, 'notransp')
%         for i1 = 1:s
%         row_ind = (i1-1)*nb+1:i1*nb;
%         tempz = ifft(fft(flip(v(row_ind,:)), n) .* fft(M(i1,:).', n));
%         z =tempz(nb:nd)+z;
%         end
        V = reshape(v,[nb s]);
        Temp = ifft(fft(flip(V),n).*fft(M.',n));
        z = sum(Temp(nb:nd,:),2);
    else
%         z = zeros(nb*s,1);
%         for i1 = 1:s
%         row_ind = (i1-1)*nb+1:i1*nb;
%         tempz = ifft(fft(flip(v), n) .* fft(M(i1,:)', n));
%         z(row_ind,1) = tempz(na:nd);
%         end
        Z = ifft(bsxfun(@times,fft(M',n),fft(flip(v),n))) ; %n*s
        z = vec(Z(na:nd,:)); %nb*s
    end
end

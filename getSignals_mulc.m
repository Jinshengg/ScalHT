function [fs, cs, H, X_star, Xs,obs] = getSignals_mulc(r, s, n,m, sep)
% input s : channel number
%       r : model order, level of spectral sparsity
%       n : length 
%       m:  number of elements  
%       sep: seperation 
% output fs: frequenc location
%         cs: amp withrespect to freq
%         H :channel coefficient
%         X_star: true siganl
%         Xs :observed siganl
%         obs::obs matrix
%-------------------------------------------------------------------
% fs = rand(r,1);
fs = rand(r,1)-0.5;
if sep    
if r == 1
    sep = 1;
else
    fs_sort = sort(fs);
    sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
end

while(sep<16/n) %primal 1/n
    fs = rand(r,1);
    if r == 1
        sep = 1;
    else
        fs_sort = sort(fs);
        sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
    end
end

end

dynamic_range= 20;
cs = exp(-1i*2*pi*rand(r,1)) .* (1 + 10.^(rand(r,1).*(dynamic_range/20)));

% dynamic_range= 10;

Cs = diag(cs); %J x J diag matrix

% channel coefficient
H = zeros(s, r);
for kk = 1:r
    H(:,kk) = randn(s,1);
    H(:,kk) = H(:, kk)/norm(H(:,kk));
end

A = zeros(r, n);
for kk  = 1:r
    tmp = amplitude(fs(kk), n);
    A(kk,:) = tmp.';
end
X_star = H * Cs * A;

K=randsample(n*s,m);
ind = zeros (n*s,1);
ind (K) = 1;
obs = reshape (ind,[s,n]);
% get the observation


Xs = X_star.*obs;

end

%%
function a = amplitude(fs, n)
N = 0:1:n-1;
a = exp(1i*2*pi*fs.*N);
end

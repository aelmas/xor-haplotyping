function [Perror, swr, err_xor] = BitPerform(Inference, TrueSolution, worstcase, Lmiss)
% [Perror, SwitchRate, err_xor] = BitPerform(Inference, TrueSolution, worstcase, Lmiss)
X = mod(sum(TrueSolution,3),2);
S = size(TrueSolution,1);
N = size(TrueSolution,2);
if nargin > 3
    missingSNPs = sum(sum(Lmiss==1));
else missingSNPs = 0;
end
if nargin < 3;
    worstcase = sum(floor(S*N-missingSNPs/2));
end
misinferred = 0;
missed = [];
swr = [];
err_xor = [];
for i = 1:N
    heterozygous    = X(:,i)==1; % heterozygous sites
    homozygous      = X(:,i)==0; % XOR sites (hidden homozygous)
    switches     = min(sum( [...
        (Inference(:,i,2) ~= TrueSolution(:,i,1))  ...
        (Inference(:,i,2) ~= TrueSolution(:,i,2))] ,1));
    if switches > 0
        missed = [missed i]; misinferred = misinferred + 1;
    end
    if sum(heterozygous) > 0
        switches_hetero = min(sum( [...
            (Inference(heterozygous,i,2) ~= TrueSolution(heterozygous,i,1))  ...
            (Inference(heterozygous,i,2) ~= TrueSolution(heterozygous,i,2))] ,1));
        swr =     [swr switches_hetero/(sum(heterozygous)/2)];
    end
    if sum(homozygous) > 0
        switches_xor = min(sum( [...
            (Inference(homozygous,i,2) ~= TrueSolution(homozygous,i,1))  ...
            (Inference(homozygous,i,2) ~= TrueSolution(homozygous,i,2))] ,1));
        err_xor = [err_xor switches_xor/sum(homozygous)];
    end % XOR sites
end
Perror  =   misinferred / N;
swr     =   mean(swr);
err_xor =   mean(err_xor);

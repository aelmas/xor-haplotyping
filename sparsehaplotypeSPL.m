function [solution, ZD, i_, PLtime] = sparsehaplotypeSPL(dataG, W, interest,Lmiss,flagBlock)
% [Inference, Dictionary, blockpartitions, PLtime] = sparsehaplotypeSPL(dataG, W, interest,Lmiss,flagBlock)
% dataG = [Genotype,Haplotype1,Haplotype2] %dims: (L)x(N)x(3)
% L = max number of SNPs
% N = number of individuals
% W = max partition size
if nargin < 3; interest = 1:size(dataG,2); end
if nargin < 2; L = size(dataG,1); W = L; end
solution = [];
flag_fetch = 0;                                               
N = size(dataG,2);
L = size(dataG,1);
tic;
if flagBlock==1
    display('PL method')
    [i_, sub_partitions, ZD_partitions] = ...
        blockPartition(dataG(:,:,1),W,Lmiss(:,interest));       
elseif flagBlock==2
    display('fixed block partitioning')
    flag_fetch = 0;
    i_ = [];
    for j=1:ceil(L/W)
        i_ = [i_; (W*(j-1)+1) min(W*j,L)];
    end
elseif flagBlock==0    
    i_ = [1 L]
end
PLtime = toc
% PL METHOD                                                   
Q = size(i_,1);
[lP,lB] = listProcess(Q);
n_ligations = length(lP);
for Ligate = 1: n_ligations
    solution = [];
    Dictionary = [];
    cols = [];
    rows = [];
    b_size = 2^(Ligate-1);
    if Ligate == n_ligations; b_size = Q; end
    n_pairs = lP(Ligate);
    n_blocks = lB(Ligate);
    for s = 1:n_blocks
        block = i_((s-1)*b_size+1,1):i_(min(s*b_size,Q),2)
        % Partition Ligation
        sub_data_int = dataG(block,:,1);
        sub_Lmiss    = Lmiss(block,interest);
        if Ligate == 1
            Restricted_Z = getZ_xor(sub_data_int);
            if flag_fetch == 1
                % fetch first inference for efficiency
                sub_solution = sub_partitions(block,:,:);
                ZD = ZD_partitions(block,:);
            else
                % recompute inference
                sub_data_int(sub_Lmiss==1) = 6;
                [estf,H12] = fSubMod(sub_data_int,Restricted_Z);
                sub_solution = zeros(size(sub_data_int,1),size(sub_data_int,2),2);
                sub_solution(:,:,1) = Restricted_Z(:,H12(:,1));
                sub_solution(:,:,2) = Restricted_Z(:,H12(:,2));
                ZD = uniondata(sub_solution);
            end
        else
            ss = (s)*2;
            Restricted_Z = concatenate( R_Z(:,:,ss-1:min(ss,lB(Ligate-1))),...
                M_C(ss-1:min(ss,lB(Ligate-1))), M_R(ss-1:min(ss,lB(Ligate-1))) );
                sub_data_int(sub_Lmiss==1) = 6;
            [estf,H12] = fSubMod(sub_data_int,Restricted_Z);
            sub_solution = zeros(size(sub_data_int,1),size(sub_data_int,2),2);
            sub_solution(:,:,1) = Restricted_Z(:,H12(:,1));
            sub_solution(:,:,2) = Restricted_Z(:,H12(:,2));
            ZD = uniondata(sub_solution);
        end
        solution = [solution ;sub_solution];
        Dictionary(1:size(ZD,1),1:size(ZD,2),s) = ZD;
        cols(s) = size(ZD,2);
        rows(s) = size(ZD,1);
    end %Each block
    M_C = cols;
    M_R = rows;
    R_Z = Dictionary;
end %Ligation Phases

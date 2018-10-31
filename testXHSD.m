clear all; close all;
% % individuals     = str2num( argv(){1} )
% % nFullGpercent   = str2num( argv(){2} )
% % per_miss        = str2num( argv(){3} )
% % per_errtyped    = str2num( argv(){4} )
% % W               = str2num( argv(){5} )
% % max_realizations= str2num( argv(){6} )
%-------------------------------------------------------------------------
% % Parameters
%-------------------------------------------------------------------------
individuals         = [100]   % number of individuals in the CFTR dataset
nFullGpercent       = [50]    % percentage of the individuals full-genotyped
per_miss            = [0]     % percentage of the missing data
per_errtyped        = [0]     % percentage of typing errors
W                   = 8     % max block size
max_realizations    = 10    % number of repeated runs to average out noise
flagBlock = 1;  % 1:PL method, 2:fixed-length partitioning, 0:no-partition
flagError = 1; 
%-------------------------------------------------------------------------
% % CFTR data
%-------------------------------------------------------------------------
dataHchar =    ['00001110001101010010000'; ...
                '00001110001101010010111'; ...
                '10101110001101010010000'; ...
                '00010000001101010010000'; ...
                '00010110001101010010000'; ...
                '00001010001101010010000'; ...
                '10101000001101010010000'; ...
                '11010001001101010010000'; ...
                '00010110001101010010111'; ...
                '11001000001101010010111'; ...
                '00010000001101010010100'; ...
                '00011101001101010010100'; ...
                '00001000001101010010111'; ...
                '10101000001101010010111'; ...
                '10101000110111010011111'; ...
                '00010001110010101000111'; ...
                '00101000110010101000000'; ...
                '00001001000010101000111'; ...
                '10010110110010101000101'; ...
                '00001000110010101001100'; ...
                '10101001001101001110111'; ...
                '10101111001101001110111'; ...
                '10101001001101001110000'; ...
                '11001001001101001110000'; ...
                '10110000001101001110111'; ...
                '10110000001101001110000'; ...
                '10001001001101001100000'; ...
                '10101001110010001000000'; ...
                '10101000110111001000000']';
% % % dataHchar = ['000';'001';'010';'000';'001';'010']'
dataH = [];
for i = 1:size(dataHchar,2)
    haplotype_i = ...
        logical(str2num(flipud([fliplr(dataHchar(:,i)');...
        blanks(length(dataHchar(:,i)'))]'))')';
    dataH = [dataH haplotype_i];
    
end
[dataGuniq,dataGiden] = gen_xor_genotypes(dataH);
dataG = [dataGuniq dataGiden]; dataX = dataG(:,:,1);
% %-----------------------------------------------------------------------
meanPerror = [];
meanSwitchRate = [];
meanError_xor = [];
cardinality_D = [];
nFullG = round(individuals*nFullGpercent/100);
Perror = []; 
SwitchRate = [];
Error_xor = [];
HeterozygousRate = [];
cardinality_vec = [];
runtimes = [];
noMTI = 0;
MTI_indexes = []; divlist = [];
for realization = 1:max_realizations
    solution = []; ZD = [];
    % Missing data
    Pmiss = per_miss/100;
    Perrtyped = per_errtyped/100;
    [Xmiss, Missing_entries] = pick_data(dataX,Pmiss);
    [~, typing_errors] = pick_data(dataX,Perrtyped);
    L = size(dataG,1); N = size(dataG,2);
    % MTI
    Xs = []; MTI_indexes = []; qt = 0;
    while length(MTI_indexes) < nFullG && qt<1000        
        qt = qt+1;
        indices = pickrandomX(1:size(dataG,2),size(dataG,2));
        [xs, mti_indexes] = findMTI(dataG(:,indices,1));
        Xs = uniondata([Xs xs]);
        MTI_indexes = uniondata([MTI_indexes indices(mti_indexes)]);
    end
    if length(MTI_indexes)<nFullG
        disp(['WARNING: # full genotypes found: ' ...
            num2str(length(MTI_indexes)) ', needed: ' num2str(nFullG)])        
    else
        disp(['# informative individuals found: ' num2str(length(MTI_indexes))]);
        (MTI_indexes)
    end
    interest = [MTI_indexes ...
        pickrandomX(setdiff(1:size(dataG,2),MTI_indexes),...
        individuals-length(MTI_indexes))];
    X           = dataX(:,interest,1);      %shuffle signals
    gr_solution = dataG(:,interest,2:3);    %shuffle signals
    Lmiss = Missing_entries;
    Lmiss_interest = Lmiss(:,interest);
    typing_errors_interest = typing_errors(:,interest);    
    Xsumsort = sum(sort(Xs,2),1);
    if nFullGpercent~=0 && Xsumsort(1)~=0
        disp('Error Xs: no subset of genotypes can resolve all the SNPs');
    end
    if isempty(MTI_indexes); noMTI = noMTI+1;
    end
    G = sum(gr_solution(:,1:length(MTI_indexes),1:2),3);
    %----------------------------------------------
    % Augment with full genotypes
    %----------------------------------------------
    %             0: 00     homozygous-common 0
    %             2: 11     homozygous-mutant 1
    %             1: 01/10  heterozygous
    %             4: 00/11  homozygous-hidden (XOR-site)
    %             6: missing SNP
    XG = X;
    if flagError
        XG(typing_errors_interest==1) = ...
            randi(2,size(XG(typing_errors_interest==1)))-1;
    end
    XG(XG==0) = 4;
    XG(:,1:length(MTI_indexes)) = G;
    %----------------------------------------------
    % Sparse Haplotyping with PL Method
    %----------------------------------------------
    tic;
    display('starting XHSD...');
    [solution, ZD, partitions, PLtime] = ...
        sparsehaplotypeSPL(XG, W, interest, Lmiss, flagBlock);
    elapsedtime = toc;
    [Pe,SR,err_xor] = BitPerform(solution, gr_solution);
    %----------------------------------------------
    Perror     = [Perror Pe];
    SwitchRate = [SwitchRate SR];
    Error_xor  = [Error_xor err_xor];
    cardinality_vec = [cardinality_vec ...
        size(uniondata([solution(:,:,1) solution(:,:,2)]),2)];
    runtimes   = [runtimes elapsedtime];
end
meanPerror     = mean(Perror)
meanSwitchRate = mean(SwitchRate)
meanError_xor  = mean(Error_xor)
cardinality_D  = mean(cardinality_vec);
meanRuntime    = mean(runtimes)
save(['XHSD_results_' num2str(individuals) ...
    '_' num2str(nFullGpercent) '_' num2str(per_miss) '_' num2str(per_errtyped)])
%
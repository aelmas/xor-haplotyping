function [listProc,listBloc] = listProcess(Q)
listProc = Q;
listBloc = Q;
n_blocks = Q;
while n_blocks > 1
    n_merges = floor(n_blocks/2);
    n_left = mod(n_blocks,2);
    n_blocks = n_merges + n_left;
    listProc = [listProc n_merges];
    listBloc = [listBloc n_blocks];
end


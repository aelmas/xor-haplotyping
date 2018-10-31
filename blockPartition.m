function [i_out, sub_partitions, ZD_partitions] = blockPartition(Genotypes,W,Lmiss)
C = [];
i_stars = [];
filled = [];
i_out = [];
L = size(Genotypes,1);
sub_partitions  = zeros(size(Genotypes));
sub_partitions1 = zeros(size(Genotypes));
sub_partitions2 = zeros(size(Genotypes));
ZD_partitions = [];
for k = 1:L;
    display([num2str(k) '/' num2str(L) ' th total entropy is being computed']);
    [c,i_s,sub_hap] = getC_i_k(C,k,Genotypes,W,Lmiss);
    C(k) = c;
    i_stars(k) = i_s;
    sub_partitions1(i_s:k,:,k) = sub_hap(:,:,1);
    sub_partitions2(i_s:k,:,k) = sub_hap(:,:,2);
end
while L > 0
    indx = i_stars(L);
    filled(1:L) = indx;
    L = indx-1;
end
un = uniondata(filled);
for ind = un
    arr = find(filled==ind);
    i_out = [i_out ; [arr(1) arr(end)]];
    sub_partitions(arr,:,1) = sub_partitions1(arr,:,arr(end));
    sub_partitions(arr,:,2) = sub_partitions2(arr,:,arr(end));
end
ZD_partitions = uniondata(sub_partitions);
end % function blockPartition

function [c,i_s,hap_is_k] = getC_i_k(C,k,Genotypes,W,Lmiss)
statement = [];
ii = max(1,k-W+1);
[E_ii_k, hap_ii_k] = getE_i_k(ii,k,Genotypes,Lmiss);
statement(ii) = 0 + E_ii_k;
for i = ii+1:k;
    [E_i_k, sub_hap] = getE_i_k(i,k,Genotypes,Lmiss);
    statement(i) = C(i-1) + E_i_k;
end
[c,ii_s] = min(statement(ii:k));
i_s = ii_s + ii-1;
hap_is_k = hap_ii_k(ii_s:end,:,:);
end % function getC_i_k

function [E_i_k,Haplotypes] = getE_i_k(i,k,Genotypes,Lmiss)
block = Genotypes(i:k, :); Lmiss_block = Lmiss(i:k, :);
%
Restricted_Z = getZ_xor(block);
block(Lmiss_block==1) = 6;
[estf,H12] = fSubMod(block,Restricted_Z);
Haplotypes = zeros(size(block,1),size(block,2),2);
Haplotypes(:,:,1) = Restricted_Z(:,H12(:,1));
Haplotypes(:,:,2) = Restricted_Z(:,H12(:,2));
f = estf(estf~=0);
E_i_k = -(f'*log2(f));
end % function getE_i_k

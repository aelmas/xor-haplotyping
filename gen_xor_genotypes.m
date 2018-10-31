function [dataG,recomb] = gen_xor_genotypes(dataH)
% [dataG, recomb] = gen_xor_genotypes(dataH)
% [[.. gi ..](dim1),[.. h1i ..](dim2),[.. h2i ..](dim3)] (dims=LxNx3) = gen_genotypes(dataH)
dataG = [];
recomb = [];
[L,M] = size(dataH);
i = 0;
k = 0;
for l = 1:M
    for m = l:M
        h1i = dataH(:,l);
        h2i = dataH(:,m);
        gi = mod(h1i+h2i,2);
        same = 0;
        for j = 1:i
            if isequal(gi,dataG(:,j,1))
                same = 1; break;
            end
        end
        if same == 0
            i = i+1;
            dataG(:,i,1) = gi;
            dataG(:,i,2) = h1i;
            dataG(:,i,3) = h2i;
        else
            k = k+1;
            recomb(:,k,1) = gi;
            recomb(:,k,2) = h1i;
            recomb(:,k,3) = h2i;
        end
    end
end


function Z = getZ_xor(dataset)
%Z = getZ_xor(dataset)
[L,N] = size(dataset);
Z = [];
for i = 0:2^L-1
    Z = [Z dec2binvec(i,L)'];
end

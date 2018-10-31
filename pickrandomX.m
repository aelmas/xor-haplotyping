function picked = pickrandomX(vector,X)
% picked = pickrandomX(vector,X)
lV = length(vector);
% if X >= lV; display('Warning: All the genotypes are picked'); end
pick = randperm(lV);
picked = vector(pick(1:min(X,lV)));

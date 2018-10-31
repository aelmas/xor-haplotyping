function [X_pick, Lpick, loci_pick]= pick_data(X, Ppick)
[L,N] = size(X);
X_pick = {};
loci_pick = {};
Lpick = zeros(L,N);
for i = 1:N
    gi = [];
    for j = 1:L
        if rand >= Ppick %don't pick
            gi = [gi; X(j,i)];
        else
            loci_pick = [loci_pick; [i j]];
            Lpick(j,i) = 1;
        end
    end
    X_pick = [X_pick gi];
end


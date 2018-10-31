function U = uniondata(u,option)
% U = uniondata(u,option)
% u input
% option = 1, remove same columns
% option = 2, remove same rows
if nargin < 2; option = 1;
elseif option == 2; u = u';
end
[row,col] = size(u);
if col == 0
    U = u;
else U = u(:,1);
end
if col > 1;
    for i = 2:col
        same = false;
        for j = 1:size(U,2)
            if sum(u(:,i) == U(:,j)) == row
                same = true; break
            end
        end
        if ~same
            U = [U u(:,i)];
        end
    end
end
if option == 2; U = U'; end

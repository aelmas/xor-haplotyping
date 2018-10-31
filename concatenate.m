function conted = concatenate(data,cols,rows)
conted = [];
c = 0;
if size(data,3) == 1
    conted = data(1:rows(1),1:cols(1));
else
    primary = data(1:rows(1),:,1);
    secondary = data(1:rows(2),:,2);
    for i = 1:cols(1)
        for j = 1:cols(2)
            c = c + 1;
            conted(:,c) = [ primary(:,i) ; secondary(:,j) ];
        end
    end
end
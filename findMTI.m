function [Xs, indexes] = findMTI(X)
% function [Xs, indexes] = findMTI(X) 
% O(n k) time complexity
% n cols, # of genotypes
% k rows, # of SNPS
Xs = [];
indexes = [];
X = X';
[n,k] = size(X);
%n rows are genotyptes
%k columns are SNPs
other_Xi = 1:n;
% 1. Apply Theroem 2
% a. choose a spine X1
for spine = randi(n);
% spine = 1;
mask_other_Xi = ones(1,n);
mask_other_Xi(spine) = 0;
other_Xi = other_Xi(mask_other_Xi==1);
X1 = X(spine,:);
nonspine_edge_labels = (1:k) .* ~X1;
% b.
for i = other_Xi
    Xi = X(i,:);
    if sum(X1.*Xi) == 0
        Xs = [X1; Xi]';
        indexes = [spine; i]';
        return 
    end
end
% c.
% d.
footprints = [];
intervals  = [];
for i = other_Xi
    Xi = X(i,:);
    footprints(i,:) = Xi.*X1;
    X1_n_Xi_labels = (1:k) .* (Xi.*X1);
    li = min(X1_n_Xi_labels);
    ri = max(X1_n_Xi_labels);
    intervals(i,:) = [li ri]; % [begins ends]
end
% e.
[r, sigma] = min(intervals(:,2)); %leftmost ending edge (SNP)
[l, ro   ] = max(intervals(:,1)); %rightmost beginning edge
% f.
if r < l-1
    Xs = [X(sigma,:); X(ro,:)]';
    indexes = [sigma; ro]';
    return
end
% g.
if r == l-1; % v = [r l]; % pivot vertex
% 2. Analyze the pivot vertex: 1 ... (o) r  (v)  l (o)
% a. If there is a path Xi ending at v:
for i = other_Xi
    % If Xi includes l, (begins from the rightmost beginning, l),
    % return {Xsigma Xi}
    if  intervals(i,1) == l
        Xs = [X(sigma,:); X(i,:)]';
        indexes = [sigma; i]';
        return
        % otherwise,        (ends at the leftmost end, r),
        % return {Xi Xro}
    elseif  intervals(i,2) == r
        Xs = [X(i,:); X(ro,:)]';
        indexes = [i; ro]';
        return
    end
end
% % b. If degree(v) = 3, (Fig. 7c),
%
% %% LOOK FOR EVERY Xi 
% % Xi include r, Xj include l,
% % Xi Xj intersect before v,
% % return {X1 Xsigma Xro)
% for i = other_Xi
%     % Xi includes r but l
%     if footprints(i,r) == 1 && footprints(i,l) == 0
%         for j = other_Xi;
%             if j ~= i;
%                 if footprints(j,l) == 1 && footprints(j,r) == 0
%                     
%                 end
%             end
%         end
%     end
% end
%
% %% OR check the tails of Xsigma and Xro
% %% if one path's tail covers the other's
if footprints(sigma,r) == 1 % && 
    tail_size = min(length(r+1:n),length(1:l-1));
    if isequal(X(sigma,r+1:r+1+(tail_size-1)), X(ro,l-1-(tail_size-1):l-1))
        Xs = [X1; X(sigma,:); X(ro,:)]';
        indexes = [spine; sigma; ro]';
        return
    end 
end
% c. Otherwise, let e be a nonspine edge incident on v, (Fig. 7d):
% pick e randomly; e = nonspine_edge_labels(randi(length(nonspine_edge_labels))); 
    % i.
    for e = nonspine_edge_labels
        for i = other_Xi
            for j = other_Xi
            if  j ~= i
                if X(i,l) == 1 && X(i,e) == 1
                if X(j,l) == 0 && X(j,e) == 0
                    Xs = [X(i,:); X(j,:)]';
                    indexes = [i; j]';
                    return
                end                
                end
            end
            end
        end
    end    
    % ii.
    for e = nonspine_edge_labels
        for i = other_Xi
            for j = other_Xi
            if  j ~= i
                if X(i,r) == 1 && X(i,e) == 1
                if X(j,r) == 0 && X(j,e) == 0
                    Xs = [X(i,:); X(j,:)]';
                    indexes = [i; j]';
                    return
                end                
                end
            end
            end
        end
    end    
    
end %pivot analysis       
end %spine
return

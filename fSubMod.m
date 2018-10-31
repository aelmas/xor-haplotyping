function [estf,haplotype2] = fSubMod(g, phi, NumIterations)
% g is base-4 matrix:
%
% 0: 00     homozygous-common-0 (G)
% 2: 11     homozygous-mutant-1 (G)
% 1: 01/10  heterozygous        (G or X)
% 4: 00/11  homozygous-hidden   (XOR-site) (X)
% 6: missing
bi=4;
NumberPeople = size(g,2);
haplotype = zeros(NumberPeople,2);
haplotype2 = zeros(NumberPeople,2);
for person = 1:NumberPeople
    Ambiguities(person)=length(union(find(g(:,person)==1),find(g(:,person)==4))); 
    Missing(person)=length(find(g(:,person)==6));
end
D = [];
Columns1=find(Ambiguities==0);
Columns2=find(Missing);
Columns = setdiff(Columns1,Columns2);
for i=1:length(Columns)
    for j = 1:size(phi,2)
        if(isequal(g(:,Columns(i))/2,phi(:,j)))
            if(ismember(j,D))
                haplotype(Columns(i),1)=find(D==j);
                haplotype(Columns(i),2)=find(D==j);
            else
                D = [D j];
                haplotype(Columns(i),1)=length(D);
                haplotype(Columns(i),2)=length(D);
            end
        end
    end
end
NoD = 1:size(phi,2); NoD = setdiff(NoD, D);
ToResolve = union(find(Ambiguities),find(Missing));
Resolved = intersect(find(Ambiguities==0),find(Missing==0));
% Compute L0
for i=1:NumberPeople
    indexes{i} = setdiff(1:length(g(:,i)),find(g(:,i)==6));
    hidd{i}      = find(g(:,i)==4);
    xi = g(:,i);
    xi(hidd{i}) = 0;
    xi = xi(indexes{i});
    Ls(i)=norm(xi)^2;
end
L0 = mean(Ls);
flag =1;
counter=1;
maxValue(1)=0;
AddedPhi = [];
Indicator = [];
for i=1:length(D)
    for j=i:length(D)
        AddedPhi = [AddedPhi phi(:,D(i))+phi(:,D(j))];
        Indicator = [Indicator;i j];
    end
end
NumberPeople = size(g,2);
% Compute Projections
minError = zeros(NumberPeople,1);
tempToResolve = ToResolve;
for person=1:length(ToResolve)
    flag = 1;
    gi = g(:,ToResolve(person));
    gi(hidd{ToResolve(person)}) = 0;
    gi = gi(indexes{ToResolve(person)});
    weightXg = (bi*norm(double(~(mod(g(:,ToResolve(person)),2)).*g(:,ToResolve(person))~=4),1)/length(indexes{ToResolve(person)}));
    for i=1:size(AddedPhi,2)
        ri = AddedPhi(:,i);
        ri(hidd{ToResolve(person)}) = mod(ri(hidd{ToResolve(person)}),2);
        ri = ri(indexes{ToResolve(person)});
        temp(i)=weightXg*norm(gi-ri,2)^2;                                               
        %pdist([g(:,ToResolve(person))';AddedPhi(:,i)'],'chebychev')^2;%norm(g(:,ToResolve(person))-AddedPhi(:,i))^2;
        if(temp(i)==0)
            if(flag)
                flag=0;
                tempToResolve = setdiff(tempToResolve, ToResolve(person));
                haplotype(ToResolve(person),1)=Indicator(i,1);
                haplotype(ToResolve(person),2)=Indicator(i,2);
            end
        end
    end
    if(isempty(D))
        xi = g(:,ToResolve(person)); xi(hidd{ToResolve(person)})=0;
        minError(ToResolve(person))=weightXg*norm(xi)^2;                                
    else
        [minValue, indicator(ToResolve(person))]=min(temp);
        minError(ToResolve(person))=minValue;
    end
end
ToResolve = tempToResolve;
flag = 1;
if(isempty(ToResolve))
    flag = 0;
    Testf = zeros(length(D),1);
    for person=1:NumberPeople
        Testf((haplotype(person,1)))=Testf((haplotype(person,1)))+1;
        Testf((haplotype(person,2)))=Testf((haplotype(person,2)))+1;
        haplotype2(person,1)=D(haplotype(person,1));
        haplotype2(person,2)=D(haplotype(person,2));
    end
    Testf= Testf/sum(Testf);
    estf = zeros(size(phi,2),1);
    estf(D)=Testf;
end
while(flag)
    resultVsNoD = zeros(length(NoD),1);
    for i=1:length(NoD)
        Matrix = [];
        Indicator = [];
        for j=1:length(D)
            Matrix = [Matrix phi(:,D(j))+phi(:,NoD(i))];
        end
        Matrix = [Matrix phi(:,NoD(i))+phi(:,NoD(i))];
        result = zeros(NumberPeople,1);
        for person=1:length(ToResolve)
            temp=[];
            gi = g(:,ToResolve(person));
            gi(hidd{ToResolve(person)}) = 0;
            gi = gi(indexes{ToResolve(person)});
            weightXg = (bi*norm(double(~(mod(g(:,ToResolve(person)),2)).*g(:,ToResolve(person))~=4),1)/length(indexes{ToResolve(person)}));
            wX = (length(indexes{ToResolve(person)})^2)/length(g(:,ToResolve(person)))^2;
            for j=1:size(Matrix,2)
                Mi = Matrix(:,j);
                Mi(hidd{ToResolve(person)}) = mod(Mi(hidd{ToResolve(person)}),2);
                Mi = Mi(indexes{ToResolve(person)});
                temp(j)=weightXg*wX*norm(gi-Mi)^2;    
            end
            [minValue, indicator]=min(temp);
            result(ToResolve(person))=minValue;
            if(result(ToResolve(person))<minError(ToResolve(person)))
            else
                result(ToResolve(person))= minError(ToResolve(person));
            end
        end
        resultVsNoD(i) = mean(result); 
    end
    [maxValueTemp, posMax] = max(L0 - resultVsNoD);
    possibilities = find((L0 - resultVsNoD)==maxValueTemp);
    if (length(possibilities)~=1)
        posMax = possibilities(randi(length(possibilities)));
    end
    AddedPhi = [];
    Indicator = [];
    for j=1:length(D)
        AddedPhi = [AddedPhi phi(:,D(j))+phi(:,NoD(posMax))];
        Indicator = [Indicator;j length(D)+1];
    end
    AddedPhi = [AddedPhi phi(:,NoD(posMax))+phi(:,NoD(posMax))];
    Indicator = [Indicator;length(D)+1 length(D)+1];    
    tempToResolve = ToResolve;
    for person=1:length(ToResolve)
        flag2 = 1;
        temp=[];
        gi = g(:,ToResolve(person));
        gi(hidd{ToResolve(person)}) = 0;
        gi = gi(indexes{ToResolve(person)});
        weightXg = (bi*norm(double(~(mod(g(:,ToResolve(person)),2)).*g(:,ToResolve(person))~=4),1)/length(indexes{ToResolve(person)}));
        for i=1:size(AddedPhi,2)
            ri = AddedPhi(:,i);
            ri(hidd{ToResolve(person)}) = mod(ri(hidd{ToResolve(person)}),2);
            ri = ri(indexes{ToResolve(person)});
            temp(i)=weightXg*norm(gi-ri)^2;                                             
            if(temp(i)==0)
                if(flag2)
                    flag2=0;
                    tempToResolve = setdiff(tempToResolve, ToResolve(person));
                    haplotype(ToResolve(person),1)=Indicator(i,1);
                    haplotype(ToResolve(person),2)=Indicator(i,2);
                end
            end
        end
        [minValue, indicator]=min(temp);
        result(ToResolve(person))=minValue;
        if(result(ToResolve(person))<minError(ToResolve(person)))
            minError(ToResolve(person))=result(ToResolve(person));
        else
            result(ToResolve(person)) = minError(ToResolve(person));
        end
    end
    ToResolve = tempToResolve;    
    maxValue(counter)=maxValueTemp;
    D = [D NoD(posMax)];
    NoD = setdiff(NoD, NoD(posMax));
    counter = counter + 1;
    if(isempty(ToResolve))
        flag = 0;
        Testf = zeros(length(D),1);
        for person=1:NumberPeople
            Testf((haplotype(person,1)))=Testf((haplotype(person,1)))+1;
            Testf((haplotype(person,2)))=Testf((haplotype(person,2)))+1;
            haplotype2(person,1)=D(haplotype(person,1));
            haplotype2(person,2)=D(haplotype(person,2));
        end
        Testf= Testf/sum(Testf);
        estf = zeros(size(phi,2),1);
        estf(D)=Testf;
    end
end

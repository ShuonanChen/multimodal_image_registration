function [trans_set,trans_set_lower_bound,max_viol,max_trans_set]=maximal_translation_group(X,Y,margin);
% Code to find the maximal translation set
% Code not double checked yet
DX=squareform(pdist(X));
DY=squareform(pdist(Y));
x_idx=cell(size(X,1),size(Y,1));
y_idx=cell(size(X,1),size(Y,1));
for i=1:size(X,1)
    for j=1:size(Y,1)
        [II0,JJ0]=find(abs(DX(i,:)'-DY(j,:))<margin);
        keep=find((II0~=i).*(JJ0~=j));
        x_idx{i,j}=II0(keep);
        y_idx{i,j}=JJ0(keep);
    end
end
trans_set=zeros(size(X,1),size(Y,1));
trans_set_lower_bound=zeros(size(X,1),size(Y,1));
max_viol=zeros(size(X,1),size(Y,1));
max_trans_set=cell(size(X,1),size(Y,1));
tic;t=0;numsearch=size(X,1)*size(Y,1);
for i=1:size(X,1)
    for j=1:size(Y,1)
        t=t+1;
        if length(x_idx{i,j})>trans_set_lower_bound(i,j)
            T=Y(j,:)-X(i,:);
            finalI=x_idx{i,j};finalJ=y_idx{i,j};%all pairs of points in X,Y whose distances to x_i and y_j are epsilon close
            tmp=(norms(X(finalI,:)+T-Y(finalJ,:),2)<margin);
            finalI=finalI(tmp);finalJ=finalJ(tmp);
            x_set=[i finalI']';
            y_set=[j finalJ']';
            ind=sub2ind(size(trans_set),x_set,y_set);
            trans_set_lower_bound(ind)=max(trans_set_lower_bound(ind),length(x_set));
            trans_set(ind)=max(trans_set(ind),trans_set_lower_bound(ind));
        end
        clc
        fprintf(['Searching for transformations (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc;
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
        
    end
end
end

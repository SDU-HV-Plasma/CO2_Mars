function [ B,G,Y,Bsim,Gsim ] = HeO_matrix( X,nsp,number_of_reactions,Ralldlt,sim,Vibration,n,number_of_reaction )

global B G Bsim Gsim

% 建立有向关系
% 建立相互关系矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        B{i,j}=intersect(find(X(i,1:number_of_reactions)>0),find(X(j,1:number_of_reactions)<0));
    end
    i=i+1;
end

B{3,1}=[6,8,9,12,B{3,1}];
B{4,1}=[12,B{4,1}];
B{6,1}=[7,B{6,1}];
B{9,1}=[1,B{9,1}];
B{11,1}=[7,8,13,B{11,1}];
B{3,2}=[34,B{3,2}];
B{4,2}=[16,34,B{4,2}];
B{5,2}=[17,B{5,2}];
B{7,2}=[15,B{7,2}];
B{6,4}=[14,B{6,4}];

if Vibration==1
 
    Bgg=cell(length(n),length(n));
    Bgv=cell(length(n),nsp-length(n));
    Bvg=cell(nsp-length(n),length(n));
    Bvv=cell(nsp-length(n),nsp-length(n));
    Bgg{2,1}=[2,4];
    Bgg{3,1}=[8,12];
    Bgg{4,1}=9;
    Bgg{9,1}=3;
    Bgg{10,1}=5;
    Bgg{11,1}=11;
    Bgg{4,2}=23;
    Bgg{3,4}=16;
    Bvg{1,1}=1;
    Bvg{3,1}=[8,13];
    Bvg{4,1}=10;
    Bvg{4,3}=26;
    
    for i=1:length(n)
        for j=1:length(n)
            Bgg{i,j}=Bgg{i,j}+number_of_reaction;
            B{i,j}=union(B{i,j},Bgg{i,j},'legacy');
        end
    end
    
     for i=1:length(n)
        for j=1:nsp-length(n)
            Bgv{i,j}=Bgv{i,j}+number_of_reaction;
            B{i,j+length(n)}=union(B{i,j+length(n)},Bgv{i,j},'legacy');
        end
     end
    
     for i=1:nsp-length(n)
        for j=1:length(n)
            Bvg{i,j}=Bvg{i,j}+number_of_reaction;
            B{i+length(n),j}=union(B{i+length(n),j},Bvg{i,j},'legacy');
        end
     end
    
     for i=1:nsp-length(n)
        for j=1:nsp-length(n)
            Bvv{i,j}=Bvv{i,j}+number_of_reaction;
            B{i+length(n),j+length(n)}=union(B{i+length(n),j+length(n)},Bvv{i,j},'legacy');
        end
     end
    
end

% 拓扑矩阵G
G=B;%针对三体进行补充
G{2,1}=[15,G{2,1}];
G{4,1}=[14,G{4,1}];
G{2,3}=[16,17,G{2,3}];
G{4,3}=[14,19,G{4,3}];
G{1,4}=[9,G{1,4}];
G{2,4}=[15,17,G{2,4}];
G{3,4}=[18,G{3,4}];
G{1,5}=[12,G{1,5}];
G{2,5}=[34,G{2,5}];
G{2,9}=[35,G{2,9}];

if Vibration==1
 
    Ggg=cell(length(n),length(n));
    Ggv=cell(length(n),nsp-length(n));
    Gvg=cell(nsp-length(n),length(n));
    Gvv=cell(nsp-length(n),nsp-length(n));
    Ggg{1,2}=1;
    Ggg{1,3}=13;
    Ggg{1,4}=[8,10];
    Ggv{1,1}=2;
    Ggv{2,1}=29;
    Ggv{1,2}=4;
    Ggv{1,3}=[12,16];
    Ggv{1,4}=9;
    Ggv{2,4}=23;
    Ggv{1,5}=26;
    
    for i=1:length(n)
        for j=1:length(n)
            Ggg{i,j}=Ggg{i,j}+number_of_reaction;
            G{i,j}=union(G{i,j},Ggg{i,j},'legacy');
        end
    end
    
     for i=1:length(n)
        for j=1:nsp-length(n)
            Ggv{i,j}=Ggv{i,j}+number_of_reaction;
            G{i,j+length(n)}=union(G{i,j+length(n)},Ggv{i,j},'legacy');
        end
     end
    
     for i=1:nsp-length(n)
        for j=1:length(n)
            Gvg{i,j}=Gvg{i,j}+number_of_reaction;
            G{i+length(n),j}=union(G{i+length(n),j},Gvg{i,j},'legacy');
        end
     end
    
     for i=1:nsp-length(n)
        for j=1:nsp-length(n)
            Gvv{i,j}=Gvv{i,j}+number_of_reaction;
            G{i+length(n),j+length(n)}=union(G{i+length(n),j+length(n)},Gvv{i,j},'legacy');
        end
     end
    
end

for i=1:nsp
    for j=1:nsp            
        B{i,j}=setdiff(B{i,j},Ralldlt);
        G{i,j}=setdiff(G{i,j},Ralldlt);
    end
end

% if Vibration~=1
%     for i=1:nsp
%         for j=1:nsp
%             if i>length(n)+1 || j>length(n)+1
%                 B{i,j}=[];
%                 G{i,j}=[];
%             end
%         end
%     end
% 
%     for i=1:length(n)
%         for j=1:length(n)
%             B{i,j}(B{i,j}>number_of_reaction)=[];
%             G{i,j}(G{i,j}>number_of_reaction)=[];
%         end
%     end
% end

% 建立有向矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        Y1(i,j)=sum(X(i,B{i,j}));
    end
    Y2(i)=sum(X(i,X(i,:)>0));
    i=i+1;
end
Y=Y1./Y2';

Bsim=B;Gsim=G;

if sim==1
    Bsim(find(all(cellfun(@(x) isempty(x),Bsim),2)),:)=[];
    Bsim(:,find(all(cellfun(@(x) isempty(x),Bsim))))=[];
    Gsim(find(all(cellfun(@(x) isempty(x),Gsim),2)),:)=[];
    Gsim(:,find(all(cellfun(@(x) isempty(x),Gsim))))=[];
%     Y(find(all(cellfun(@(x) isnan(x),Y),2)),:)=[];
%     Y(:,find(all(cellfun(@(x) isempty(x),Y))))=[];
end

end


function [fastestPath_Distance,fastestPath_type,fastestPath_Rate,Distance,Predistance,Impt_Species] = Dosimplification(Ratesumexini,nsp,species,sim,nall,closeness,betweenness,dt,B,number_of_reactions,X,Rall)
global Ratesumexini nsp species sim nall closeness betweenness dt B number_of_reactions X Rall 
[fastestPath_Distance, fastestPath_type, Distance, Predistance]=fastestPath(Ratesumexini, nsp, species); %��������,�����ǳ�ʼ�ڵ�Ũ��
if sim==1
    nall(nall==1e-3)=[];
    fastestPath_type(find(all(cellfun(@(x) isempty(x),fastestPath_type),2)),:)=[];
    fastestPath_type(:,find(all(cellfun(@(x) isempty(x),fastestPath_type))))=[];
    fastestPath_Distance(logical(eye(size(fastestPath_Distance))))=inf;
    fastestPath_Distance(all(isinf(fastestPath_Distance),2),:)=[];
    fastestPath_Distance(:,all(isinf(fastestPath_Distance),1))=[];
end
fastestPath_Rate = log10(nall)-fastestPath_Distance ;%a.u.%��������
fastestPath_Rate(fastestPath_Rate==-Inf)=0;
fastestPath_Rate=fastestPath_Rate-min(min(fastestPath_Rate));
fastestPath_Rate(logical(eye(size(fastestPath_Rate))))=0;
Impt_Species = fastestPath_Rate * (closeness'* 0.01 + betweenness * 0.2 ) .* log10(nall) ;

skeleton_reaction=[];
m1=['�������豣����������Ŀ:'];
m1=input(m1);
[~,psp]=sort(Impt_Species,'descend');
M1=cellstr(species(psp(1:m1)));
disp(['����Ҫ�������Ϊ:',M1]);
r(1:number_of_reactions)=1;
Bs=B;
%���漰����ɾȥ���ӵķ�Ӧɾ�������γ�ʣ���������ܹ��ɵ�����Bs
for i=1:nsp
    for j=1:nsp
        if ismember(i,sort(psp(m1+1:nsp))) || ismember(j,sort(psp(m1+1:nsp)))
            r(B{i,j})=0;
%             Bs{i,j}=[];
        end
    end
end
reaction_delete=find(r==0);
for i=1:nsp
    for j=1:nsp
        Bs{i,j}=setdiff(Bs{i,j},reaction_delete);
    end
end
max_reaction=sum(r==1);
n_max_reaction=find(r==1);
for i=1:nsp
    if ismember(i,sort(psp(1:m1))) && (all(cellfun(@isempty,Bs(:,i))==1) || all(cellfun(@isempty,Bs(i,:))==1) )
        error('����������������٣�����Ⱥ�޷���ɿɿ��ıջ���');
    end
end

m2=['�������豣���ķ�Ӧ��Ŀ:'];
m2=input(m2);
if m2>max_reaction
    error('��Ҫ��ķ�Ӧ��Ŀ���ڸ�����Ⱥ�ܹ���ɵ����Ӧ��Ŀ��');
end
% ���Bs�������������Ӽ��������ķ�Ӧ���ͼ�������
Ratemaxs=zeros(nsp,nsp);
Ratemaxs_reaction=cell(nsp,nsp);
for i=1:nsp
    for j=1:nsp
        if ~isempty(max(Rall(Bs{i,j})))
            Ratemaxs(i,j)=max(X(i,Bs{i,j}));
            Ratemaxs_reaction{i,j}=find(X(i,:)==Ratemaxs(i,j));
        end
    end
end
% Ѱ�һ�ԭ�Ǽ�
n_skeleton_reaction=[];
Ratemaxexinis=Ratemaxs*1/10*dt./nall';
Ratemaxmax1s=Ratemaxexinis;Ratemaxmax2s=Ratemaxexinis;
for i=1:nsp
    for j=1:nsp
        Ratemaxmax1s(i,Ratemaxmax1s(i,:)<max(Ratemaxmax1s(i,:)))=0;
        Ratemaxmax2s(Ratemaxmax2s(:,j)<max(Ratemaxmax2s(:,j)),j)=0;
    end
end
skeleton=Ratemaxmax1s;
skeleton(Ratemaxmax1s~=Ratemaxmax2s)=Ratemaxmax1s(Ratemaxmax1s~=Ratemaxmax2s)+Ratemaxmax2s(Ratemaxmax1s~=Ratemaxmax2s);
skeleton_reaction=Ratemaxs_reaction;
for i=1:nsp
    for j=1:nsp
        if skeleton(i,j)==0
            skeleton_reaction{i,j}=[];
        end
        n_skeleton_reaction=[n_skeleton_reaction skeleton_reaction{i,j}];
    end
end
n_skeleton_reaction=unique(n_skeleton_reaction);
skeleton_num=length(n_skeleton_reaction);
if m2<skeleton_num
    error('��Ҫ��ķ�Ӧ��Ŀ�������Ļ�ԭ�Ǽܣ�')
end
%��Bs��ȥ����ԭ�Ǽ��еķ�Ӧ���γ�Bss��Ѱ����Ҫ��ķ�Ӧ����
Bss=Bs;
Rss=cell(nsp,nsp);
n_Bss=[];n_Rss=[];reaction_sim=[];
reaction_to_remain=m2-skeleton_num;
for i=1:nsp
    for j=1:nsp
        Bss{i,j}=setdiff(Bss{i,j},skeleton_reaction{i,j});
        Rss{i,j}=Rall(Bss{i,j})./nall(j);
        n_Bss=[n_Bss Bss{i,j}];
        n_Rss=[n_Rss Rss{i,j}];
    end
end
lookup_reactions(1,:)=n_Bss;%����ʣ�෴Ӧ���ɵĲ��Ҿ��󣬵�һ��Ϊ��Ӧ���࣬�ڶ���Ϊ��Ӧ��Լ������
lookup_reactions(2,:)=n_Rss;
lookup_reactions=sortrows(lookup_reactions',2,'descend')';
for k=1:length(lookup_reactions)
    reaction_sim=[n_skeleton_reaction lookup_reactions(1,1:k)];
    if length(unique(reaction_sim))==m2
        break
    end
end
disp(['��Ӧ�ļ򻯷�Ӧ����Ϊ�� ',num2str(unique(reaction_sim)),' ,���з�Ӧ ',num2str(n_skeleton_reaction),' �����˼򻯼��Ļ�ԭ�Ǽ�']);

end

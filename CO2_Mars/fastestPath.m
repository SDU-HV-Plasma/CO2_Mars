function [fastestPath_Distance, fastestPath_Type, Distance, Predistance] = fastestPath(Ratesumexini, nsp, species)
global Ratesumexini nsp species 
Predistance=1./Ratesumexini';
Distance=log10(Predistance);
% Distance=Distance-min(min(Distance))+1;
[i,j,v] = find(Distance);
sparse_matrix = sparse(i, j, v, nsp, nsp);%构建稀疏矩阵
fastestPath_Distance = zeros(nsp, nsp);%最速路径长度
fastestPath_Type = {[]};
fastestPath_Type = repmat(fastestPath_Type, nsp, nsp);%最速路径类型
for m=1:nsp
    for n=1:nsp 
        if n~=m
            [fastestPath_Distance(m,n), fastestPath_Type{m,n}, ~] = graphshortestpath(sparse_matrix, m, n);
%             fastestPath_Distance(m,n) = fastestPath_distance;
%             fastestPath_Type{m,n} = fastestPath_type;
        end
    end
end

% bg=biograph(sparse_matrix,species,'ShowWeights','on','NodeAutoSize','off','ArrowSize',10,'EdgeFontSize',15);
% 
% m1=['请输入起始反应物：'];
% m1=input(m1);
% M1=find(strcmp(species,m1));
% m2=['请输入终末生成物：'];
% m2=input(m2);
% M2=find(strcmp(species,m2));
% [a,b,c]=graphshortestpath(sparse_matrix, M1, M2);
% 
% dolayout(bg);
% 
% for i=1:length(bg.Edges) %修改边的属性
%     bg.Edges(i).LineColor=[0 0.4470 0.7410];
%     bg.Edges(i).LineWidth=1;
% end
% 
% set(bg.Nodes(b),'Color',[1 0.4 0.4]);
% edges = getedgesbynodeid(bg,get(bg.Nodes(b),'ID'));
% set(edges,'LineColor',[1 0 0]);
% set(edges,'LineWidth',2);
% 
% nx=zeros(1,length(bg.nodes));
% ny=zeros(1,length(bg.nodes));
% for i=1:length(bg.nodes) %修改节点属性
%     bg.Nodes(i).Size=[200 60];
%     bg.Nodes(i).LineColor=[1 1 1];
%     bg.Nodes(i).Color=[0.9290 0.6940 0.1250];
%     bg.Nodes(i).FontSize=50;
%     nx(i)=bg.Nodes(i).Position(1);%储存节点的横坐标信息
%     ny(i)=bg.Nodes(i).Position(2);
% end
% 
% x=randperm(36);
% x=reshape(x,[6,6]);
% % for i=[2 3 4] %修改同级节点相对位置
% %     temp=sort(nx(1,j:j+i-1));
% %     for k=1:length(temp)
% %         bg.Nodes(j+k-1).Position(1)=temp(k);
% %     end
% %     j=j+i;
% % end
% for i=1:6
%     bg.Nodes(x(i,1)).Position(1)=2000;
%     bg.Nodes(x(i,2)).Position(1)=3100;
%     bg.Nodes(x(i,3)).Position(1)=4200;
%     bg.Nodes(x(i,4)).Position(1)=5300;
%     bg.Nodes(x(i,5)).Position(1)=6400;
%     bg.Nodes(x(i,6)).Position(1)=7500;
% end
% 
% for j=1:6
%     bg.Nodes(x(1,j)).Position(2)=100;
%     bg.Nodes(x(2,j)).Position(2)=600;
%     bg.Nodes(x(3,j)).Position(2)=1100;
%     bg.Nodes(x(4,j)).Position(2)=1600;
%     bg.Nodes(x(5,j)).Position(2)=2100;
%     bg.Nodes(x(6,j)).Position(2)=2600;
% end
% 
% % for i=unique(ny) %ny是一个存储了所有节点纵坐标的矩阵
% %     for j=1:length(bg.nodes)
% %         if bg.Nodes(j).Position(1,2)==i
% %             bg.Nodes(j).Position(1,2)=0.06*length(bg.nodes)*i;
% %         end
% %     end
% % end
% 
% dolayout(bg,'Pathsonly',true);
% 
% view(bg);
end


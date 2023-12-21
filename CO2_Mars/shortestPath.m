function minWeightMatrix=shortestPath(G,nodeNum)
 
minWeightMatrix=zeros(nodeNum,nodeNum);
for i=1:nodeNum %获取每个顶点到其它节点最短路径权值，Dijkstra 源点遍历
    D=G(i,:);   %初始化节点到各个邻居节点的权值
    D(D==0)=inf;
    tag=zeros(1,nodeNum); %0 为对应位置节点未求出最大路径
    tag(i)=1;
    while sum(tag,2)~=nodeNum
        noTag=find(tag==0); %还未获得最短距离的节点
        minValue=min(D(noTag));
        if minValue~=inf
            tempIndex=find(D==minValue); %找出候选节点，更新值
            tag(tempIndex)=1;
            % 将候选节点去掉
            for candidateNode=1:length(tempIndex)
                noTag(noTag==tempIndex(candidateNode))=[];
            end
            % 修改源点到各个未标记的顶点间的距离
            for j=1:length(noTag) %对所有未获得最短距离的节点进行更新
                minTemp=D(noTag(j));
                for k=1:length(tempIndex)
                    if minTemp>(minValue+G(tempIndex(k),noTag(j)))&&G(tempIndex(k),noTag(j))~=0
                        minTemp=minValue+G(tempIndex(k),noTag(j));
                    end
                end
                D(noTag(j))=minTemp;
            end            
        else
            break;
        end
    end
    minWeightMatrix(i,:)=D;
end
minWeightMatrix(isinf(minWeightMatrix))=0;
minWeightMatrix=minWeightMatrix-diag(diag(minWeightMatrix));



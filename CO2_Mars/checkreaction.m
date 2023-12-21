function [check_reaction]=checkreaction(Gsim,nsp,spc_sim)

% m = ['请输入需查看所参与反应的粒子：'];
% m = input(m);
% m = strcmp(species,m);

nsp=nsp-spc_sim;

for k=1:nsp
    M1 = [];
    M2 = [];
    for n=1:nsp
        M1 = [M1,Gsim{k,n}];
%     end
%     for n=1:nsp
        M2 = [M2,Gsim{n,k}];
    end
    check_reaction{k} = unique([M1,M2]);
end
% Reaction = num2str(reaction);
% disp(strcat('所参与的的反应为：',Reaction));

end


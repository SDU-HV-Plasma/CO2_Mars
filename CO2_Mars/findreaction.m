function findreaction(B,species)

global B species

m1=['请输入第一个反应物：'];
m1=input(m1);
M1=strcmp(species,m1);
N1=unique([B{:,M1}]);
m2=['请输入第二个反应物：'];
m2=input(m2);
M2=strcmp(species,m2);
N2=unique([B{:,M2}]);
m3=['请输入第三个反应物：'];
m3=input(m3);
Thirdspecies=strcmp(m3,'null');
if(~Thirdspecies)
    M3=strcmp(species,m3);
    N3=unique([B{:,M3}]);
    N12=intersect(N1,N2);
    Reaction=intersect(N12,N3);
else
    Reaction=intersect(N1,N2);
end
Reaction=num2str(Reaction);

if(~isempty(Reaction))
    disp(strcat('符合要求的反应为：',Reaction));
else
    disp('无法找到符合要求的反应！');
end

end


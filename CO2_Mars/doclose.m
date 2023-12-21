function [ clo,num_of_clo ] = doclose( Gsim,nsp,spc_sim )

nsp=nsp-spc_sim;
clo=zeros(nsp);
num_of_clo=zeros(nsp);
% clo、num_of_clo的灵敏度均只到4的的最短路径
i=1;
while i<nsp+1
    for j=i+1:nsp
        if ~isempty(Gsim{i,j}) || ~isempty(Gsim{j,i})
            clo(i,j)=1;
            num_of_clo(i,j)=1;
        end
    end
    i=i+1;
end

clo=clo+clo';
num_of_clo = num_of_clo + num_of_clo';

i=1;
while i<nsp+1
    for j=i+1:nsp
        if clo(i,j)==0
            jump=0;
            for k=1:nsp
                if clo(k,i)==1 && clo(k,j)==1
                    jump=jump+1;
                end
            end
            if jump~=0
                clo(i,j)=2;
                num_of_clo(i,j)=jump;
            end
        end
    end
    i=i+1;
end
    
clo = triu(clo,0) + tril(clo',-1);
num_of_clo = triu(num_of_clo,0) + tril(num_of_clo',-1);

i=1;
while i<nsp+1
    for j=i+1:nsp
        if clo(i,j)==0
            jump=0;
            for k=1:nsp
                if  clo(k,i)==1 && clo(k,j)==2  
                    jump = jump + num_of_clo(k,j);
                end
            end
            if jump~=0
                clo(i,j)=3;
                num_of_clo(i,j)=jump;
            end
        end
    end
    i=i+1;
end
    
clo = triu(clo,0) + tril(clo',-1);
num_of_clo = triu(num_of_clo,0) + tril(num_of_clo',-1);

i=1;
while i<nsp+1
    for j=i+1:nsp
        if clo(i,j)==0
            jump=0;
            for k=1:nsp
                if  clo(k,i)==1 && clo(k,j)==3
                    jump = jump + num_of_clo(k,j);
                end
            end
            if jump~=0
                clo(i,j)=4;
                num_of_clo(i,j)=jump;
            end
        end
    end
    i=i+1;
end
    
clo = triu(clo,0) + tril(clo',-1);
clo(clo==0) = 5;
clo = clo - diag(diag(clo));

num_of_clo = triu(num_of_clo,0) + tril(num_of_clo',-1);
num_of_clo(num_of_clo==0) = 99;%表示需要人工帮助
num_of_clo = num_of_clo - diag(diag(num_of_clo));

end




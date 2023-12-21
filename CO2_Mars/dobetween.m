function [ Bet ] = dobetween( clo,nsp,num_of_clo,spc_sim )

nsp=nsp-spc_sim;
bet=zeros(nsp,nsp);
for m=1:nsp
    i=1;
    while i<nsp+1
        for j=i:nsp
            if clo(i,j)==2
                p=find(clo(i,m)==1 && clo(m,j)==1);
                bet(i,j)=size(p,1);
            end
            if clo(i,j)==3
                p1=0;
                p2=0;
                if clo(i,m)==2 && clo(m,j)==1
                    p1 = num_of_clo(i,m);
                end
                if clo(i,m)==1 && clo(m,j)==2
                    p2 = num_of_clo(m,j);
                end
                bet(i,j) = p1 + p2 ;
            end
            if clo(i,j)==4
                p1=0;
                p2=0;
                p3=0;
                if clo(i,m)==3 && clo(m,j)==1
                    p1 = num_of_clo(i,m);
                end
                if clo(i,m)==2 && clo(m,j)==2
                    p2 = num_of_clo(i,m).* num_of_clo(m,j);
                end
                if clo(i,m)==1 && clo(m,j)==3
                    p3 = num_of_clo(m,j);
                end
                bet(i,j) = p1 + p2 + p3;
            end
        end
        i=i+1;
    end

Bet{m,1}=bet;

end


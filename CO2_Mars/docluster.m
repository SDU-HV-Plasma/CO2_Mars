function [ clustering ] = docluster( nsp,clo,spc_sim )

nsp=nsp-spc_sim;
clustering = zeros(nsp,1);

m=1;
while m<nsp+1 
    clo_species=find(clo(m,:)==1);
    gi=size(clo_species,2);
    g=gi.*(gi-1);
    %bi=sum(sum(clo(clo_species,clo_species)))./2;
    bi=size(find(clo(clo_species,clo_species)==1),1);
    clustering(m,1)=bi./g;
    
    m=m+1;
end

end
function [dne]=electron_part_balance(Te,P,n,area_volume,Tg,rates,e_isp)
    [nsp,species,X]=rates(Te,P,n,area_volume,Tg,1);

    if(Te<Tg)
        Te=Tg;
    end

    dne=abs(sum(X(e_isp,:))); 


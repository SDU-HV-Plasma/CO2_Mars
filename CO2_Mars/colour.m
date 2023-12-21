function [cc]=colour(i)



colors=['r' 'g' 'b' 'm' 'c' 'k' 'y' ':' '--' '-.k' '-.b' '-.m' '-.y'];
N=length(colors);
k=rem(i,N);
if(~k) 
    k=N; 
end
cc=colors(k);

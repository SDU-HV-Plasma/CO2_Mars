function [K]=lookups(filename,Te)
%[K]=lookup(filename,Te)
%Returns rate coefficent K for a given Te
%Requires look up table file

aux=load(['C:\zthmatlabhome\MATLAB2\' filename]);
temp=aux(:,1);                          %create array of teperatures.
rate=aux(:,2);                          %create array of rates.

if(Te<min(temp))
    K=rate(1);
elseif (Te>max(temp))
    K=rate(end);    
else
    K=interp1(temp,rate,Te);
end
% Te(Te<min(temp))=min(temp);
% Te(Te>max(temp))=max(temp);
% K=interp1(temp,rate,Te);
end


function [ti]=tiempo(s,t)
ti=1;
for i=1:length(t)
    if s>t(i)
       ti=ti+1;
    else
    end
end
if ti>length(t)
    ti=length(t);
end
if ti<1
    ti=1;
end

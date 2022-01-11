function snnew = mutation(snew,pmutation)
Bitlength = size(snew,2);
snnew = snew;
pmm = IfCroIfMut(pmutation);
if pmm ==1
    chb = round(rand*(Bitlength-1))+1;  
    snnew(chb) = abs(snew(chb)-1);
end
end
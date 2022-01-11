function scro = crossover(population,seln,pc)
Bitlength = size(population,2);
pcc = IfCroIfMut(pc);  
if pcc ==1
    chb = round(rand*(Bitlength-2))+1; 
    scro(1,:) = [population(seln(1),1:chb) ...
        population(seln(2),chb+1:Bitlength)];
    scro(2,:) = [population(seln(2),1:chb)...
        population(seln(1),chb+1:Bitlength)];
else
    scro(1,:) = population(seln(1),:);
    scro(2,:) = population(seln(2),:); 
end
end
function[Fitvalue,cumsump] = fitnessfun(layer,population,boundsbegin,...
    boundsend,Bitlength_conductivity,Bitlength_thick,SimColle)
popsize = size(population,1); 
population_conductivity = cell(1,layer);    
population_thick = cell(1,layer-1);         
Fitvalue = zeros(1,popsize);               
for j = 1: popsize
    for i=0:layer-1
        population_conductivity{i+1} = population(:,...
            i*Bitlength_conductivity+1:...
            (i+1)*Bitlength_conductivity);
        if i<layer-1
            population_thick{i+1} = population(:,...
                (layer*Bitlength_conductivity+1+i*Bitlength_thick):...
                layer*Bitlength_conductivity+(i+1)*Bitlength_thick);
        end
        x.conductivity(j) = transform2to10(population_conductivity{i+1}(j,:));
        if i<layer-1
            x.thick(j) = transform2to10(population_thick{i+1}(j,:));
        end 
        xx.conductivity(i+1) = boundsbegin.conductivity+x.conductivity(j)*...
            (boundsend.conductivity-boundsbegin.conductivity)/...
            (power(2,Bitlength_conductivity)-1);
        if i<layer-1
            xx.thick(i+1) = boundsbegin.thick+x.thick(j)*...
                (boundsend.thick-boundsbegin.thick)/...
                (power(2,Bitlength_thick)-1);
        end        
    end
    [IP,QP]= Forward_FDEM(xx.conductivity,xx.thick);
    Fitvalue(j) = targetfun(IP,QP,SimColle);
end
Fitvalue = Fitvalue';
%% Fitness value sensitivity reconstruction
X_Index = var(Fitvalue);
if X_Index>1e12
    Index =1; 
elseif X_Index<=1e12 && X_Index>1e10
    Index =0.8;
elseif X_Index<=1e10 && X_Index>1e6
    Index =0.6;
elseif X_Index<=1e6 && X_Index>1e3
    Index =0.4;
elseif X_Index<=1e3 && X_Index>1e1
    Index =0.3;
else
    Index =0.2;
end
fsum = sum(Fitvalue);
Pperpopulation = Fitvalue/fsum;
Pperpopulation = (1-Pperpopulation).^1; % (Index*popsize)
Pfsum = sum(Pperpopulation);
Pperpopulation = Pperpopulation/Pfsum;
cumsump(1) = Pperpopulation(1);
%-------------------------------------------------------------------------%
for j = 2:popsize
    cumsump(j) = cumsump(j-1)+Pperpopulation(j);
end
cumsump = cumsump';
end

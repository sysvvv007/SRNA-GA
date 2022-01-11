%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application of sensitivity reconstruction nonlinear adaptive improved
% genetic algorithm for frequency domain electromagnetic inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%% Initialization
%
clc;  close all; clear;
dbstop if error
addpath(genpath('functions'))
global freq

%
%% Imput
%
%---------------------------- Emission frequency -------------------------%
LowestFreq = 90;                                                           % Lowest frequency£¨Hz£©
HighestFreq = 96e3;                                                        % Highest frequency£¨Hz£©
%---------------------------- Geophysical model --------------------------%
Real.conductivity =  [0.01,5,0.01];                                        % Conductivity (S/m)
Real.thick =   [5,6];                                                      % Thickness of layers (m)
bounds.conductivity = [0 20];                                              % Range of inversion conductivity
bounds.thick = [0.1 20];                                                   % Range of layers thickness
SimError = 0;                                                              % Noise level
N.Inversion = 5;                                                             % Inversion times 50
N.population = 20;                                                         % Multiple of first population 20

%
%% Correlation calculation
%
syms XfreqSmall XfreqBig
freqLow = 10^XfreqSmall- LowestFreq;
freqHigh = 10^XfreqBig- HighestFreq;
[freqLow] = solve(freqLow,XfreqSmall);
freqLow = eval(freqLow);
[freqHigh] = solve(freqHigh,XfreqBig);
freqHigh = eval(freqHigh);
freq = logspace(freqLow,freqHigh,10);
layer = length(Real.conductivity);
[Real.IP,Real.QP] = Forward_FDEM(Real.conductivity,Real.thick);
SimColle.IP = Real.IP+SimError.*(-1+2.*rand(size(Real.IP))).*Real.IP;
SimColle.QP = Real.QP+SimError.*(-1+2.*rand(size(Real.QP))).*Real.QP;

%
%% Parameters of improved genetic algorithm
%
Generationmax = 300;
pcrossover1 = 0.99;
pcrossover2 = 0.6;
pmutation1 = 0.5; 
pmutation2 = 0.2; 
precision = 0.001;
popsize =ceil(0.0005*2*layer*((bounds.thick(:,2)-...
    bounds.thick(:,1)))/precision);

%
%% Improved genetic algorithm
%
boundsbegin.conductivity = bounds.conductivity(:,1);
boundsend.conductivity = bounds.conductivity(:,2);
boundsbegin.thick = bounds.thick(:,1);
boundsend.thick = bounds.thick(:,2);
%-------------------------------------------------------------------------%
Bitlength_conductivity = ceil(log2((boundsend.conductivity-...
    boundsbegin.conductivity)'./precision));
Bitlength_thick = ceil(log2((boundsend.thick-...
    boundsbegin.thick)'./precision));
Bitlength = layer*Bitlength_conductivity+(layer-1)*Bitlength_thick;
%-------------------------------------------------------------------------%
gpuArray.Con1 = zeros(1,N.Inversion);
gpuArray.Con2 = zeros(1,N.Inversion);
gpuArray.Con3 = zeros(1,N.Inversion);
gpuArray.thk1 = zeros(1,N.Inversion);
gpuArray.thk2 = zeros(1,N.Inversion);
gpuArray.ConAsm = cell(1,N.Inversion);
gpuArray.thkAsm = cell(1,N.Inversion);
gpuArray.IP_Mul = cell(1,N.Inversion);
gpuArray.QP_Mul = cell(1,N.Inversion);
gpuArray.MulEvl = zeros(1,N.Inversion);
gpuArray.time = zeros(1,N.Inversion);
gpuArray.MulEvl_Ground = zeros(1,N.Inversion);
%-------------------------------------------------------------------------%
working=waitbar(0,'Please wait....');
for D=1:N.Inversion
%% Bayesian statistical optimization scheme for first population 
    tic;
    population = round(rand(N.population*popsize,Bitlength));
    [Fitvalue,~] = fitnessfun(layer,population,boundsbegin,...
        boundsend,Bitlength_conductivity,Bitlength_thick,SimColle);
    Fitvalue_Matrix = [Fitvalue';1:N.population*popsize];
    SortFit = sortrows(Fitvalue_Matrix',1,'ascend') ;
    BeterSortFit = SortFit(1:popsize,:);
    BeterSortFit = sortrows(BeterSortFit,2,'ascend');
    Position_FitBetter = BeterSortFit(:,2)';
    population =  population(Position_FitBetter,:);
    [Fitvalue,cumsump] = fitnessfun(layer,population,boundsbegin,...
        boundsend,Bitlength_conductivity,Bitlength_thick,SimColle);
    time_1thPopulation = toc;
%-------------------------------------------------------------------------%
    Generation = 1;
    ymin = zeros(1,Generationmax);
    ymean = zeros(1,Generationmax);
    tic;
    while Generation<Generationmax+1
        scnew = [];
        smnew = [];
        for jj = 1:2:popsize
%% Calculation of nonlinear adaptive mutation rate and crossover rate 
            seln = selection(cumsump);
            Fitvalue_AGA(1) = Fitvalue(seln(1));
            Fitvalue_AGA(2) = Fitvalue(seln(2));
            Meanfit = mean(Fitvalue);
            Smallerfit = min(Fitvalue_AGA);
            Biggerfit = max(Fitvalue_AGA);
            Minfit = min(Fitvalue);
            Maxfit = max(Fitvalue);
            if Biggerfit <= Meanfit
                pcrossover =  (Meanfit*pcrossover2 - Minfit*pcrossover1)...
                    /(Meanfit - Minfit) + (Biggerfit*(pcrossover1 -...
                    pcrossover2))/(Meanfit - Minfit);
                pmutation = pmutation1*exp(((Meanfit*...
                    log(pmutation2/pmutation1))-(Smallerfit*...
                    log(pmutation2/pmutation1)))/(Meanfit - Minfit));
            else
                pcrossover = pcrossover1;
                pmutation = pmutation1;
            end
%-------------------------------------------------------------------------%
            scro = crossover(population,seln,pcrossover);
            scnew(jj,:) = scro(1,:);
            scnew(jj+1,:) = scro(2,:);
            smnew(jj,:) = mutation(scnew(jj,:),pmutation);
            smnew(jj+1,:) = mutation(scnew(jj+1,:),pmutation);
            if Fitvalue(seln(1))== Minfit||Fitvalue(seln(2))== Minfit
                Fitvalue_seln = cat(2,Fitvalue(seln(1)),Fitvalue(seln(2)));
                PosSeln = find(Fitvalue_seln == min(Fitvalue_seln));
                if PosSeln(1) ==1
                    smnew(jj,:) = population(seln(1),:);
                else
                    smnew(jj+1,:) = population(seln(2),:);
                end
            end
        end
        population = smnew;
        [Fitvalue,cumsump] = fitnessfun(layer,population,boundsbegin,...
            boundsend,Bitlength_conductivity,Bitlength_thick,SimColle);
        [fmin,nmin] = min(Fitvalue);
        fmean = mean(Fitvalue);
        ymin(Generation) = fmin;
        ymean(Generation) = fmean;
        population_conductivity = cell(1,layer);
        population_thick = cell(1,layer-1);
        conductivity = zeros(1,layer);
        thick = zeros(1,layer-1);
        for ii=0:layer-1
            population_conductivity{ii+1} = population(:,...
                ii*Bitlength_conductivity+1:...
                (ii+1)*Bitlength_conductivity);
            if ii<layer-1
                population_thick{ii+1} = population(:,...
                    (layer*Bitlength_conductivity+1+ii*Bitlength_thick):...
                    layer*Bitlength_conductivity+(ii+1)*Bitlength_thick);
            end
            x.conductivity = transform2to10...
                (population_conductivity{ii+1}(nmin,:));
            if ii<layer-1
                x.thick = transform2to10(population_thick{ii+1}(nmin,:));
            end
            conductivity(ii+1) = boundsbegin.conductivity+x.conductivity*...
                (boundsend.conductivity-boundsbegin.conductivity)/...
                (power(2,Bitlength_conductivity)-1);
            if ii<layer-1
                thick(ii+1) = boundsbegin.thick+x.thick*...
                    (boundsend.thick-boundsbegin.thick)/...
                    (power(2,Bitlength_thick)-1);
            end
        end
        time(D)=toc;
        waitbar((((D-1)+Generation/Generationmax))/N.Inversion,working);
        Generation = Generation+1;
    end
    MulSol.conductivity{D} = conductivity;
    MulSol.thick{D} = thick;
    Con1(D) = MulSol.conductivity{D}(1);
    Con2(D) = MulSol.conductivity{D}(2);
    Con3(D) = MulSol.conductivity{D}(3);
    thk1(D) = MulSol.thick{D}(1);
    thk2(D) = MulSol.thick{D}(2);
%-------------------------------------------------------------------------%
    MulEvl_Ground(D) = abs(Con1(D)-Real.conductivity(1))+...
        abs(Con2(D)-Real.conductivity(2))+abs(Con3(D)-Real.conductivity...
        (3))+abs(thk1(D)-Real.thick(1))+abs(thk2(D)-Real.thick(2));
    ConAsm{D} = [Con1(D),Con2(D),Con3(D)];
    thkAsm{D} = [thk1(D),thk2(D)];
    [IP_Mul{D},QP_Mul{D}]= Forward_FDEM(ConAsm{D},thkAsm{D});
    MulEvl(D) = sum(abs(IP_Mul{D}-SimColle.IP)./abs(Real.IP))+...
        sum(abs(QP_Mul{D}-SimColle.QP)./abs(Real.QP));
end
close(working);

%
%% Output
%
Evl_AllGA = 1e2*N.Inversion/sum(MulEvl);
fprintf('Algorithm evaluation parameters:%10.5f\n',Evl_AllGA);            % The higher the value, the better the effect
%-------------------------------------------------------------------------%
EffTim_Ratio = 1e2*N.Inversion*Evl_AllGA/(time_1thPopulation+sum(time));
fprintf('Time efficiency evaluation parameters:%10.5f\n',EffTim_Ratio);   % The higher the value, the better the effect
%-------------------------------------------------------------------------%
figure;
subplot(1,2,1);
plot(1:length(Con1),ones(1,length(Con1)).*Real.conductivity(1),...
    "r:","lineWidth",2.2);
hold on;
grid on;
scatter(1:length(Con1),Con1,20,'o','filled','MarkerFaceColor',...
    [0 .6797 0.9336],'MarkerEdgeColor',[.1757 .4688  .8359],"lineWidth",1);
xlabel("Inversion times",'FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',12,'FontName','Times New Roman');
ylabel("Conductivity (s/m)",'FontName','Times New Roman');
set(get(gca,'YLabel'),'FontSize',12,'FontName','Times New Roman');
title("Conductivity",'FontName','Times New Roman');
xlim([0 N.Inversion]);
ylim([bounds.conductivity(1) bounds.conductivity(2)]);
%-------------------------------------------------------------------------%
subplot(1,2,2);
plot(1:length(thk1),ones(1,length(thk1)).*Real.thick(1),"r:",...
    "lineWidth",2.2);
hold on;
grid on;
scatter(1:length(thk1),thk1,20,'o','filled','MarkerFaceColor',...
    [0 .6797 0.9336],'MarkerEdgeColor',[.1757 .4688  .8359],"lineWidth",1);
xlabel("Inversion times",'FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',12,'FontName','Times New Roman');
ylabel("Thickness (m)",'FontName','Times New Roman');
set(get(gca,'YLabel'),'FontSize',12,'FontName','Times New Roman');
title("Thickness",'FontName','Times New Roman');
xlim([0 N.Inversion]);
ylim([bounds.thick(1) bounds.thick(2)]);
hold off;
%------------------ Best and worst results of inversion ------------------%
MulEvl_Nor = 1./MulEvl;
MulEvl_Nor = MulEvl_Nor./max(MulEvl_Nor);
MulEvl_Ground_Nor = 1./MulEvl_Ground;
MulEvl_Ground_Nor = MulEvl_Ground_Nor./max(MulEvl_Ground_Nor);
PosMin = find(MulEvl_Nor == min(MulEvl_Nor));
PosMax = find(MulEvl_Nor == max(MulEvl_Nor));
%-------------------------------------------------------------------------%
[Inv.IP,Inv.QP]= Forward_FDEM(MulSol.conductivity{PosMax},...
    MulSol.thick{PosMax});
Fig_IQppm(freq,[[Real.IP;Real.QP];[Inv.IP;Inv.QP]]); 
legend("Real.IP","Real.QP","Inv.IP","Inv.QP"); 
title("Best output compared with theory");
set(gca,'xscale','log');
%-------------------------------------------------------------------------%
[Inv.IP,Inv.QP]= Forward_FDEM(MulSol.conductivity{PosMin},...
    MulSol.thick{PosMin});
Fig_IQppm(freq,[[Real.IP;Real.QP];[Inv.IP;Inv.QP]]); 
legend("Real.IP","Real.QP","Inv.IP","Inv.QP"); 
title("Worst output compared with theory");
set(gca,'xscale','log');
%% Forward modeling response and earth model evaluation system 
figure;
bar([MulEvl_Ground_Nor', MulEvl_Nor']);
grid on;
legend("Ground Model","Forward Model",'FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',12);
ylabel("Evaluation value",'FontName','Times New Roman');
set(get(gca,'YLabel'),'FontSize',12);
xlabel("Inversion times",'FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',12);
title("Evaluation of Inversion",'FontName','Times New Roman');
% save Workspace
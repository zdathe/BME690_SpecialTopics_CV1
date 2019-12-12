clearvars
clc
geodata1 = getgeodata('GPL1261'); %Entire Mouse genome
geodata5 = getgeodata('GPL8217'); %42 normal human tissues including heart
geodata6 = getgeodata('GPL3738'); %10 normal dog tissues including heart

% not working, even with geoseriesread() vs getgeodata()
% geodata7 = getgeodata('GSE135805'); %Transplanted (healthy) and explanted (different diseases) human heart tissue

EHG = getgeodata('GPL570'); %Entire Human Genome

HNH_1 = getgeodata('GSM1053922'); %Human Normal heart_G16
HNH_2 = getgeodata('GSM1053929'); %Human Normal heart_G38
HNH_3 = getgeodata('GSM1053939'); %Human Normal heart_G66
HNH_4 = getgeodata('GSM1053940'); %Human Normal heart_G67
HNH_5 = getgeodata('GSM1053942'); %Human Normal heart_G75

HIC_1 = getgeodata('GSM1053914'); %Human Ischemic Cardiomyopathy_G1
HIC_2 = getgeodata('GSM1053916'); %Human Ischemic Cardiomyopathy_G5
HIC_3 = getgeodata('GSM1053920'); %Human Ischemic Cardiomyopathy_G12
HIC_4 = getgeodata('GSM1053921'); %Human Ischemic Cardiomyopathy_G15
HIC_5 = getgeodata('GSM1053923'); %Human Ischemic Cardiomyopathy_G19
HIC_6 = getgeodata('GSM1053927'); %Human Ischemic Cardiomyopathy_G30
HIC_7 = getgeodata('GSM1053928'); %Human Ischemic Cardiomyopathy_G33
HIC_8 = getgeodata('GSM1053930'); %Human Ischemic Cardiomyopathy_G42
HIC_9 = getgeodata('GSM1053931'); %Human Ischemic Cardiomyopathy_G43
HIC_10 = getgeodata('GSM1053932'); %Human Ischemic Cardiomyopathy_G53
HIC_11 = getgeodata('GSM1053934'); %Human Ischemic Cardiomyopathy_G55
HIC_12 = getgeodata('GSM1053936'); %Human Ischemic Cardiomyopathy_G60

HDC_1 = getgeodata('GSM1053915'); %Human Dilated Cardiomyopathy_G2
HDC_2 = getgeodata('GSM1053917'); %Human Dilated Cardiomyopathy_G6
HDC_3 = getgeodata('GSM1053918'); %Human Dilated Cardiomyopathy_G8
HDC_4 = getgeodata('GSM1053919'); %Human Dilated Cardiomyopathy_G9
HDC_5 = getgeodata('GSM1053924'); %Human Dilated Cardiomyopathy_G21
HDC_6 = getgeodata('GSM1053925'); %Human Dilated Cardiomyopathy_G24
HDC_7 = getgeodata('GSM1053926'); %Human Dilated Cardiomyopathy_G28
HDC_8 = getgeodata('GSM1053933'); %Human Dilated Cardiomyopathy_G54
HDC_9 = getgeodata('GSM1053935'); %Human Dilated Cardiomyopathy_G59
HDC_10 = getgeodata('GSM1053937'); %Human Dilated Cardiomyopathy_G62
HDC_11 = getgeodata('GSM1053938'); %Human Dilated Cardiomyopathy_G63
HDC_12 = getgeodata('GSM1053941'); %Human Dilated Cardiomyopathy_G68

%Train NN on HNH, HIC, and HDC, test on EHG (GWAS) to show clustering

HIC_all=[];
HDC_all=[];
for i=1:12
    a=eval(sprintf('HIC_%d',i));
    HIC_all = horzcat(HIC_all, a.Data(:,2));
    
    a=eval(sprintf('HDC_%d',i));
    HDC_all = horzcat(HDC_all, a.Data(:,2));
end

colnames = {'Var1','Var2','Var3','Var4','Var5','Var6','Var7','Var8',...
    'Var9','Var10','Var11','Var12'};
rownames = compose('%d', HDC_1.Data(:,1));
HDC_all=[array2table(HDC_all, 'RowNames', rownames,...
    'VariableNames', colnames); array2table(repmat(1, [1,12]),...
    'VariableNames', colnames)];

rownames = compose('%d', HIC_1.Data(:,1));
HIC_all=[array2table(HIC_all, 'RowNames', rownames,...
    'VariableNames', colnames); array2table(repmat(2, [1,12]),...
    'VariableNames', colnames)];

HNH_all=[];
for i=1:5
    a=eval(sprintf('HNH_%d',i));
    HNH_all = horzcat(HNH_all, a.Data(:,2));
end

colnames = {'Var1','Var2','Var3','Var4','Var5'};
rownames = compose('%d', HNH_1.Data(:,1));
HNH_all=[array2table(HNH_all, 'RowNames', rownames,...
    'VariableNames', colnames); array2table(repmat(3, [1,5]),...
    'VariableNames', colnames)];

writetable(cell2table(rownames), 'rownames.txt');
writetable(HDC_all, 'HDC_all.txt');
writetable(HIC_all, 'HIC_all.txt');
writetable(HNH_all, 'HNH_all.txt');

HIC = readtable('HIC_all.txt');
HDC = readtable('HDC_all.txt');
HNH = readtable('HNH_all.txt');
HIC=table2array(HIC);
HDC=table2array(HDC);
HNH=table2array(HNH);
rownames = compose('%d', table2array([readtable('rownames.txt');{1}]));

x=horzcat(HDC(1:length(HDC)-1,:), HIC(1:length(HIC)-1,:),...
HNH(1:length(HNH)-1,:));

% Filter out low values
[mask1, x, names] = genelowvalfilter(x, cell2mat(rownames(1:length(rownames)-1,:)), 'Percentile', 30);

% Filter out low entropy
[mask2, x, names] = geneentropyfilter(x, names,'prctile',30);

targ = horzcat(HDC(length(HDC),:), HIC(length(HIC),:),...
HNH(length(HNH),:));

% Probabilities
t = zeros(3,numel(targ));
for n = 1:numel(targ)
    t(targ(n),n) = 1;
end

genes=[];
indeces=[];
values=[];
fgenes=[];
findeces=[];
fvalues=[];
pall=[];
ptrain=[];
pval=[];
ptest=[];

fgenes = table2array(readtable('fgenes.txt'));
findeces = table2array(readtable('findeces.txt'));
fvalues = table2array(readtable('fvalues.txt'));
values=fvalues;

pall = table2array(readtable('pall.txt'));
ptrain = table2array(readtable('ptrain.txt'));
ptime = table2array(readtable('ptime.txt'));
ptest = table2array(readtable('ptest.txt'));
pval = table2array(readtable('pval.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find top 10 clusters based on MMSE
for i=1:10
    i=length(fgenes)+1;
    error=Inf;
    for j=1:length(x)
        if ~ismember(j, findeces)
            genes(i)=str2double(names(j));
            indeces(i)=j;
            values(i,:)=x(j,:);
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %NN BlackBox
            tic
            net = patternnet([5 5], 'traingd');
            [net, tr] = train(net,values,t,'useGPU','yes');
            y = net(values,'useGPU','yes');
            perf = mse(net,t,y);
            E=toc;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if perf < error   
                i+j/length(x)
                error=perf;
                
                ptime(i)=E;
                pall(i)=perf;
                fgenes(i)=genes(i);
                findeces(i)=j;
                fvalues(i,:)=values(i,:);

                % Recalculate Training, Validation and Test Performance
                trainTargets = t .* tr.trainMask{1};
                valTargets = t  .* tr.valMask{1};
                testTargets = t  .* tr.testMask{1};
                ptrain(i) = mse(net,trainTargets,y);
                pval(i) = mse(net,valTargets,y);
                ptest(i) = mse(net,testTargets,y);
            end
        end
    end
    writetable(array2table(ptime),'ptime.txt')
    writetable(array2table(pall),'pall.txt');
    writetable(array2table(fgenes),'fgenes.txt');
    writetable(array2table(findeces),'findeces.txt');
    writetable(array2table(fvalues),'fvalues.txt');
    writetable(array2table(ptrain),'ptrain.txt');
    writetable(array2table(pval),'pval.txt');
    writetable(array2table(ptest),'ptest.txt');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluate top N by running >30 times for tdist
index=50;
h1=5;
h2=5;
m=length(t(:,1));
for i=1:10
    x=fvalues(1:i,:);
    N=length(x);
    h1=round(sqrt((m+2)*N)+2*sqrt(N/(m+2)));
    h2=round(m*sqrt(N/(m+2)));
    
    for tdist=1:index
        disp(i+tdist/index)
        tic

        net = patternnet([h1 h2], 'traingd');
        [net, tr] = train(net,values,t,'useGPU','yes');
        y = net(values,'useGPU','yes');
        errors = gsubtract(t,y);
        Overall(tdist) = mse(net,t,y);

        %Percent error/accuracy
        predictor = y;
        for m=1:length(y)
            for n=1:3
                if predictor(n,m)==max(y(:,m))
                    predictor(n,m)=1;
                else
                    predictor(n,m)=0;
                end
            end
        end
        count=0;
        for m=1:length(y)
            if predictor(:,m)==t(:,m)
                count=count+1;
            end
        end
        percent(tdist) = count/length(y);
    
        % Recalculate Training, Validation and Test Performance
        trainTargets = t .* tr.trainMask{1};
        valTargets = t  .* tr.valMask{1};
        testTargets = t  .* tr.testMask{1};
        Train(tdist) = mse(net,trainTargets,y);
        Val(tdist) = mse(net,valTargets,y);
        Test(tdist) = mse(net,testTargets,y);

        ElapseTime(tdist) = toc;
    end

    perfPercent(i,:) = [mean(percent) std(percent)];
    perfOverall(i,:) = [mean(Overall) std(Overall)];
    perfTrain(i,:) = [mean(Train) std(Train)];
    perfVal(i,:) = [mean(Val) std(Val)];
    perfTest(i,:) = [mean(Test) std(Test)];
    perfTime(i,:) = [mean(ElapseTime) std(ElapseTime)];
end

writetable(array2table(perfOverall),'SFS_pall.txt');
writetable(array2table(perfTrain),'SFS_ptrain.txt');
writetable(array2table(perfVal),'SFS_pval.txt');
writetable(array2table(perfTest),'SFS_ptest.txt');
writetable(array2table(perfTime),'SFS_ptime.txt');
writetable(array2table(perfPercent), 'SFS_perc.txt');
addpath(genpath('../../'))

rng('default')
tic
levels = 0.2:0.1:0.9;

typeOneError = zeros(length(levels),1);
mcSamples=100;
numShuffles=1500;
sampleSize=1000;


pool = parpool;
parfor depLev= 1:length(levels);
    tic
    alpha=0.05;
    partialResults = zeros(mcSamples,1);
    bootstrapedValuesShift=[];
    model = arima('Constant',0,'AR',levels(depLev),'Variance',1- levels(depLev)^2);

    for i=1:mcSamples    
        processes = simulate(model,sampleSize,'numPaths',2)
        X=processes(:,1);
        Y=processes(:,2);
        sigX = median_heur(X);
        sigY = median_heur(Y);
        if mod(i-1,10)==0            
            [bootShift,bootstrapedValuesShift] = customShiftHSIC(X,Y,alpha,50,min(sampleSize,numShuffles),sigX,sigY);   
        else
            bootShift = customShiftHSIC(X,Y,alpha,50,min(sampleSize,numShuffles),sigX,sigY,bootstrapedValuesShift); 
        end       
        partialResults(i) = bootShift.areDependent; 
    end    
    toc
    typeOneError(depLev,:) = mean(partialResults,1);
end
toc
delete(pool)
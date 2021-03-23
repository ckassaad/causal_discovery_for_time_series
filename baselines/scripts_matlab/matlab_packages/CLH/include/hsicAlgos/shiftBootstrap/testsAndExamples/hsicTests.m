function tests = hsicTests
addpath('../');
tests = functiontests(localfunctions);
end

function setup(testCase)
rng('default');
ln=20;
ar = exp(-1/ln);
variance = 1-exp(-2/ln);
model = arima('Constant',0,'AR',{ar},'Variance',variance);
processes = simulate(model,500,'numPaths',3);
testCase.TestData.processes = processes;
end

function testIndependentForSimpeHSIC(testCase)
processes = testCase.TestData.processes;
test=shiftHSIC(processes(:,1),processes(:,2));
verifyEqual(testCase,test.areDependent,false);
end


function testDependentForSimpeHSIC(testCase)
processes = testCase.TestData.processes;
test=shiftHSIC(processes(:,1).*processes(:,3),processes(:,2).*processes(:,3));
verifyEqual(testCase,test.areDependent,true);
end

function testIndependentForCustomHSIC(testCase)
processes = testCase.TestData.processes;
X=processes(:,1);
Y=processes(:,2);
sig =median_heur(X);
test= customShiftHSIC(X,Y,0.06,20,490,sig,sig);
verifyEqual(testCase,test.areDependent,false);
end

function testDependentForCustomHSIC(testCase)
processes = testCase.TestData.processes;
X=processes(:,1).*processes(:,3);
Y=processes(:,2).*processes(:,3);
sig =median_heur(X);
test= customShiftHSIC(X,Y,0.06,20,490,sig,sig);
verifyEqual(testCase,test.areDependent,true);
end

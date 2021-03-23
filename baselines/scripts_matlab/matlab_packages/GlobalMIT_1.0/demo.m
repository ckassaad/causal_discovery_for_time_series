%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: demo();
%Step through demo of the GlobalMIT toolkit

function demo();

clc;

load husmeier_yeast_100;
fprintf('Data loaded...\n');

s_true=score_MIT(data,true_net);
fprintf('Score of the true network: %f \n',s_true);
createDotGraphic(true_net,nodeName,'True DBN');

fprintf('Press any key to continue...\n');
pause

alpha=0.999;
allowSelfLoop=1;
[best_net]=globalMIT(data,alpha,allowSelfLoop);
createDotGraphic(best_net,nodeName,'GlobalMIT DBN');

fprintf('Press any key to continue...\n');
pause

[best_net_exe,score,time]=globalMIT_exe(data,alpha,0);
createDotGraphic(best_net_exe,nodeName,'GlobalMIT C++ DBN');

fprintf('Press any key to continue...\n');
pause

clear clc;
load Yu_net_5;
n_state=3;

fprintf('Press any key to continue...\n');
pause

a1= myIntervalDiscretize(a1,n_state);
a2= myIntervalDiscretize(a2,n_state);
a3= myIntervalDiscretize(a3,n_state);

fprintf('Press any key to continue...\n');
pause

alpha=0.95;
tic;[best_net]=globalMIT(a1,alpha,1);toc
createDotGraphic(best_net,[],'GlobalMIT DBN');
compare_net(best_net,true_net,1)

fprintf('Press any key to continue...\n');
pause

[b,c]=multi_time_series_cat(a1,a2,a3);
alpha=0.95;
[best_net_ab]=globalMIT_ab(b,c,alpha,1);
createDotGraphic(best_net_ab,[],'GlobalMIT DBN')
compare_net(best_net_ab,true_net,1)

fprintf('Warning: Large data set. May take up to 1h. Press any key to continue...\n');
pause

alpha=0.9999;
data= myIntervalDiscretize(data,n_state);
[best_net_exe,score,time]=globalMIT_exe(data,alpha,1);
createDotGraphic(best_net_exe,[],'GlobalMIT C++ DBN');
compare_net(best_net_exe,true_net,1)

fprintf('The end...\n');


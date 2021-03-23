clear;clc
type='p';           %type of kernel ('p' for polynomial, 'g' for gaussian)
par=2;              %parameter of the model (order for polynomial kernel, or width for the gaussian kernel
N=2100;             %number of points of simulated time series
msim=1;             %order of the model in the simulations
mtot=1;
ltot=[25  500 2000];
epstot=0:0.01:1;
sigtot=[0 0.0005 0.03];
G=zeros(2,2,length(mtot),length(ltot),length(epstot),length(sigtot));G1=G;
corrtot=zeros(length(mtot),length(ltot),length(epstot),length(sigtot));
for ie=1:length(epstot)
    for is=1:length(sigtot)
        [xx,yy]=two_ccm(N,msim,1.8,0,epstot(ie),sigtot(is));
        for im=1:length(mtot)
            m=mtot(im);
            for il=1:length(ltot);
                disp([ie is im il]);
                data=[xx yy];            %data matrix must have the dimensions [n_points n_variables]
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [np nv]=size(data); %number of observations (time points) and number of variables
                n=ltot(il);             %endpoint
                data=data(1:ltot(il),:);
                G(:,:,im,il,ie,is)=causality(data,type,par,m);            %compute causality matrix
                G1(:,:,im,il,ie,is)=causality(data,type,1,m);            %compute causality matrix
                cc=corrcoef(data(:,1),data(:,2));
                corrtot(im,il,ie,is)=cc(1,2);
            end
        end
    end
end
% for i=1:length(ltot)
%     figure;
%     subplot(3,1,1);plot(epstot,squeeze(G(1,2,:,i,:,1)),'-ob');hold on;plot(epstot,squeeze(G(2,1,:,:,:,1)),'-or')
%     legend('x \rightarrow y','y \rightarrow x');
%     title(strcat('\sigma = ',num2str(sigtot(1))));
%     subplot(3,1,2);plot(epstot,squeeze(G(1,2,:,i,:,2)),'-ob');hold on;plot(epstot,squeeze(G(2,1,:,:,:,2)),'-or')
%     title(strcat('\sigma = ',num2str(sigtot(2))));
%     subplot(3,1,3);plot(epstot,squeeze(G(1,2,:,i,:,3)),'-ob');hold on;plot(epstot,squeeze(G(2,1,:,:,:,3)),'-or')
%     title(strcat('\sigma = ',num2str(sigtot(3))));xlabel('\epsilon');
% end
figure;
subplot(3,3,1);plot(epstot,squeeze(G1(1,2,:,:,:,1)),'o');title(strcat('\sigma = ',num2str(sigtot(1))));ylabel('causality 1->2, p=1')
subplot(3,3,2);plot(epstot,squeeze(G1(1,2,:,:,:,2)),'o');title(strcat('\sigma = ',num2str(sigtot(2))));
subplot(3,3,3);plot(epstot,squeeze(G1(1,2,:,:,:,3)),'o');title(strcat('\sigma = ',num2str(sigtot(3))));
subplot(3,3,4);plot(epstot,squeeze(G(1,2,:,:,:,1)),'o');title(strcat('\sigma = ',num2str(sigtot(1))));ylabel('causality 1->2, p=2')
subplot(3,3,5);plot(epstot,squeeze(G(1,2,:,:,:,2)),'o');title(strcat('\sigma = ',num2str(sigtot(2))));
subplot(3,3,6);plot(epstot,squeeze(G(1,2,:,:,:,3)),'o');title(strcat('\sigma = ',num2str(sigtot(3))));
subplot(3,3,7);plot(epstot,squeeze(corrtot(:,:,:,1)),'d');xlabel('\epsilon');ylabel('correlation')
subplot(3,3,8);plot(epstot,squeeze(corrtot(:,:,:,2)),'d');xlabel('\epsilon');
subplot(3,3,9);plot(epstot,squeeze(corrtot(:,:,:,3)),'d');xlabel('\epsilon');
legend('25 points','500 points','2000 points')
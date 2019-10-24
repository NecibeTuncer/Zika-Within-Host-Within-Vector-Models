function  total_ARE = MonteCarlo_WithinHost_SM3_1stTrimester

%Monkey393422

% clear all
% close all
% clc

numiter = 5; 


X = zeros(10,numiter);


tforward = (0:0.1:78);

tmeasure_v = [1,11,21,31,41,51,71,81,91,101,141,211,281,351,431,501,571,641,711,781];
tmeasure_n = [1,21,31,71,81,91,101,141,211,281,351,421,571,781];
     


ic = [45000000 10 0 28 1000000 0 0]; 

true_params = 1.0e+06 *[0.000000000005050   0.000000042603184   0.000019898046524...
              0.000000485060891   0.147716353104356   0.000000279383936...
              0.000000000007006   3.504265363595249   0.000000056280792...
              0.000000022048771];   


[~,y_trp] = ode23s(@(t,y)(model_Zika(t,y,true_params)),tforward,ic);


function dy = model_Zika(t,y,k)

 dy = zeros(7,1); 

 % SM3
 dy(1) = -k(1)*y(3)*y(1); %S susceptible target cells
 dy(2) = k(1)*y(3)*y(1)-k(2)*y(2)*y(4); %I infected target cells
 dy(3) = y(2) + y(7) - k(3)*y(3); %V free virus
 dy(4) = k(4)*y(2)*y(2)*y(4)/(y(2)*y(2)+k(5)) - k(6)*y(4); %N natural killer cells
 dy(5) = -k(7)*y(7)*y(5); % Number of hofbauer cells
 dy(6) = k(7)*y(7)*y(5) - k(8)*y(6);% Number of infected hofbauer cells
 dy(7) = k(9)*y(3) + y(6) - k(10)*y(7);%Virus in fetus

end


function error_in_data = err_in_data(k)


 [t,y] = ode23s(@(t,y)(model_Zika(t,y,k)),tforward,ic);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  BN = y(tmeasure_n(:),4)';
  ZV = y(tmeasure_v(:),3)'; 
  
 error_in_data = sum((BN - ydata_n).^2) + 0.00001*sum((ZV - ydata_v).^2);
 

end
 

noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

%for noisei = 1:length(noiselevel)
  
rng default
noisei = 5;
noiselev = noiselevel(noisei)

for i= 1:numiter
    i

   ytrp_v = y_trp(tmeasure_v(:),3);
   ydata_v = noiselev*randn(1,size(tmeasure_v,2)) + ytrp_v';
   
   ytrp_n = y_trp(tmeasure_n(:),4);
   ydata_n = noiselev*randn(1,size(tmeasure_n,2)) + ytrp_n';

k = true_params;

 
 lb = [0 0 0 0 0 0 1e-8 0 0 0];
 
k =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter')); 
                        
X(:,i) = k';


end


arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;

%  figure(1)
%  plot(tforward(tmeasure_v(:)),ydata_v,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%             hold on 
%  plot(tforward, y_trp(:,3), '-b','LineWidth',3)
% 
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% title ({'Free Virus In Pregnant Monkey 827577','ZIKV Infected At The 1st Trimester'},'FontWeight','normal')
% set(gca,'yscale','log','YTick',[10^1, 10^2, 10^3, 10^4, 10^5, 10^6],'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
% ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')
% 
% 
%  figure(2)
%  plot(tforward(tmeasure_n(:)),ydata_n,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%              hold on 
%   plot(tforward, y_trp(:,4), '-b','LineWidth',3)
% 
% ylabel('Natural Killer Cells','FontSize',14,'FontName','Sans-serif');
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% set(gca, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
% title ({'Natural Killer Cell Counts In Pregnant Monkey 827577'},'FontWeight','normal')
 
%end

save('MCS_SM3_1stTrimester_ARE')
end
function  total_ARE = MonteCarlo_WithinVectorModel

clear all
close all
clc

numiter = 1000; 

% conv = zeros(6,numiter);

X = zeros(6,numiter);

 days_v = [2,3,4,5,6,7,10,14];
% 
% midgutViral = [3.89665, 3.97, 5.14, 5.574977, 5.38589, 5.5399, 5.756262, 5.33324];
 days_n = [5, 6, 7, 10, 14];
% 
% salivaViral = [2.487, 4.0495, 3.5727, 6.4888, 6.785];
     
tforward = (0:0.1:14);

tmeasure_v = [21,31,41,51,61,71,101,141];
tmeasure_n = [51,61,71,101,141];
     


initial_cond = [1000 0.0]; 





true_params =  [1.26082924880034,433378.838892983,0.110713342233325,...
                19534.2728313170,1.76753893630032,6843584.00734574];



[~,y_trp] = ode23s(@(t,y)(model_within_vector(y,true_params)),tforward,initial_cond);

    ytrp_v = log10(y_trp(tmeasure_v(:),1));
    ytrp_n = log10(y_trp(tmeasure_n(:),2));


function dy = model_within_vector(y,k)
     
 dy = zeros(2,1); 

 dy(1) = k(1)*y(1)*(1-y(1)/k(2))-k(3)*y(1); % virus in the midgut
 dy(2) = k(3)*y(1)*y(1)/(k(4) + y(1)*y(1)) + k(5)*y(2)*(1-y(2)/k(6)); %virus in the salivary gland


end

function error_in_data = err_in_data(k)
 
 [t,y] = ode23s(@(t,y)(model_within_vector(y,k)),tforward,initial_cond);
                

  VM = log10(y(tmeasure_v(:),1)');
 
  VS = log10(y(tmeasure_n(:),2)');
  
 error_in_data = sum((ydata_v - VM).^2) +  sum((ydata_n - VS).^2);
end


 

noiselevel = [0, 0.01, 0.05, 0.1, 0.2];%, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:5
 
rng default
noiselev = noiselevel(noisei)

for i= 1:numiter
    i

   ydata_v = normrnd(ytrp_v',noiselev*ytrp_v');

   ydata_n = normrnd(ytrp_n',noiselev*ytrp_n');
  

k = true_params;

 
  lb = [0 0 0 0.0 0 0];
 
k =  fminsearchbnd(@err_in_data,k,lb,[]); 
                        
X(:,i) = k';

 
end



arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;



end

% figure(1)% = figure('position', [0, 0, 1200, 1400]);
%  
%  plot(days_v,midgutViral,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%              hold on 
%  plot(tforward, (log10(Y(:,1))), '-b','LineWidth',3)
% %axis([0 80 0 10^6])
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% %title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
% %set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
%  %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
% ylabel('Virus titer Log10 TCID50/ml in the midgut','FontSize',14,'FontName','Sans-serif')
% set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
%  hold off
%  
%  
%  figure(2)% = figure('position', [0, 0, 1200, 1400])
%  
%  plot(days_n,salivaViral,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%              hold on 
%  plot(tforward, (log10(Y(:,2))), '-b','LineWidth',3)
%   ylabel('Virus titer Log10 TCID50/ml in the saliva glands','FontSize',14,'FontName','Sans-serif');
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
%  hold off 
%  
%  figure(3)
%  plot(days_v,10.^midgutViral,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%              hold on
%  plot(tforward, Y(:,1), '-b','LineWidth',3)
% %axis([0 80 0 10^6])
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% %title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
% %set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
%  %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
% ylabel('Virus titer  TCID50/ml in the midgut','FontSize',14,'FontName','Sans-serif')
% set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
%   
%  figure(4)
%   plot(days_n,10.^salivaViral,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%              hold on 
%  plot(tforward, Y(:,2), '-b','LineWidth',3)
% %axis([0 80 0 10^6])
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% %title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
% %set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
%  %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
% ylabel('Virus titer  TCID50/ml in the saliva glands','FontSize',14,'FontName','Sans-serif')
% set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');

save('MCS_Within_Vector_ARE')
end
function  total_ARE = MonteCarlo_WithinHost_SM1

clear all
close all
clc

numiter = 1000; 


X = zeros(3,numiter);

tforward = (0:0.1:10);

tmeasure = [1 11 21 31 41 51 71 81 91 101];
 
%ic = [270000 1 0]; %MONKEY 393422 DATA
 
ic = [11000000 1 0]; %MONKEY 912116 DATA

%true_params = [0.000822078521600   0.319182923786664   3.167466721624296]; %MONKEY 393422 DATA fitting values
true_params = [0.000007580039113   1.519645525685333   1.519695264689408];%MONKEY 912116 DATA fitting values


[~,y_trp] = ode45(@(t,y)(model_Zika(t,y,true_params)),tforward,ic);


function dy = model_Zika(t,y,k)

 dy = zeros(3,1); 

 % Model 1
 dy(1) = -k(1)*y(3)*y(1); % S susceptible target cells
 dy(2) =  k(1)*y(3)*y(1) - k(3)*y(2); % I infected target cells
 dy(3) =  y(2) - k(2)*y(3); % V free virus


end



function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,ic);

 ZV = y(tmeasure(:),3)'; 
  
 error_in_data = sum((ZV - ydata).^2);
 

end

 

noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:5

rng default
noiselev = noiselevel(noisei);

for i= 1:numiter
    
    i

   ytrp_n = y_trp(tmeasure(:),3);
   ydata = noiselev*randn(1,size(tmeasure,2)) + ytrp_n';

k = true_params;

 
 lb = [0 0 0];
 
k =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter')); 
                        
X(:,i) = k';


end


arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;

 
%  figure(1)
%  plot(tforward(tmeasure(:)),ydata,'Marker','.','Color',[1 0 0],...
%                 'MarkerSize',30,'LineStyle','none')
%             hold on 
%  plot(tforward, y_trp(:,3), '-b','LineWidth',3)
% 
% ylabel('Viral RNA copies per milliliters of plasma','FontSize',14,'FontName','Sans-serif')
% xlabel('Days','FontSize',14,'FontName','Sans-serif' );
% title ({'Free Virus In Male Monkey 912116'},'FontWeight','normal')
% %set(gca,'YTick',[0, 5*10^5, 10*10^5,  3*10^6],'YTickLabel',...
%   %  {'0', '5x10^5', '10^6', '3x10^6'}, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
 
end
%save('MCS_SM1_ARE')
save('MCS_SM1_MONKEY_912116_ARE')

end

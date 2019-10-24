function Fitting_Zika_Within_SM1

clear all
close all
clc
format long

%MONKEY 912116 DATA
% 
% days = [0,1,2,3,4,5,6,7,8,9,10];
% 
% plasmaViral = [0,2580,228000,2520000,1340000,264000,18900,2100,2490,2440,0];
%             
% tforward = (0:0.1:10);
% 
% tmeasure = [1:10:101];
% 
% ic=[11000000 1 0];

%MONKEY 393422 DATA

% % days = [0,1,2,3,4,5,6,7,8,9,10];
% % 
% % plasmaViral = [0,12800,65800,34100,46700,47500,82800,5480,735,578,0];
% % true data
% 
days = [0,1,2,3,4,5,7,8,9,10];

plasmaViral = [0,12800,65800,34100,46700,47500,5480,735,578,0];
%data with an outlier removed
tforward = (0:0.1:10);

tmeasure = [1 11 21 31 41 51 71 81 91 101];
% 
 ic= [270000 1 0];



 


function dy = model_Zika(t,y,k)

 dy = zeros(3,1); 

 % Model 1
 dy(1) = -k(1)*y(3)*y(1); % S susceptible target cells
 dy(2) =  k(1)*y(3)*y(1) - k(3)*y(2); % I infected target cells
 dy(3) =  y(2) - k(2)*y(3); % V free virus


end



function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,ic);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  BV = y(tmeasure(:),3)'; 
  
  error_in_data = sum((BV - plasmaViral).^2);
 

end

%=========================================================================
%MONKEY 393422 DATA

% %k = [0.000176267857343   0.132387169452577  5.5]; %error 4.487*10^9
%k= [0.000300860403383   0.412550535468823   4.488296952331541];
% k=[0.000750121816215   0.350972726077666   3.418228926045940] %diff data
 %results in the following parameters
 k = [0.000822078521600   0.319182923786664   3.167466721624296]
% err=1.2*10^9  final
%==========================================================================

%MONKEY 912116 DATA
 %k =[ 0.000007546310539   1.523778176905339   1.523760140326185];
 % k =[ 0.000007580039087   1.519645515619110   1.519695275095262];

 lb = [0 0 0];

 % for i=1:5
 %     i
% [k,fval] = fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter', 'TolX', 1e-14,'TolFun',1e-14,...
%                         'MaxIter',10000,'MaxFunEvals',10000)); 
% 
%  disp(k);
%  end
 [T,Y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,ic);
 
 
 

 figure(1)
 plot(days,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward, (Y(:,3)), '-b','LineWidth',3)

ylabel('Viral RNA copies per milliliters of plasma','FontSize',14,'FontName','Sans-serif')
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Male Monkey 393422'},'FontWeight','normal')
set(gca,'YTick',[10^4,  3*10^4, 7*10^4],'YTickLabel',...
    { '10^4', '3x10^4', '7x10^4' }, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');

figure(2)
 plot(days,log10(plasmaViral),'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward, log10(Y(:,3)), '-b','LineWidth',3)
 
  figure(3)
 plot(days,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward, (Y(:,3)), '-b','LineWidth',3)

ylabel('Viral RNA copies per milliliters of plasma','FontSize',14,'FontName','Sans-serif')
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Male Monkey 912116'},'FontWeight','normal')
end

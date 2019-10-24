function Fitting_Zika_Within_SM2_Monkey393422

clear all
close all
clc
format long

%MONKEY 393422 DATA MALE MONKEY

days = [0,1,2,3,4,5,6,7,8,9,10];

plasmaViral = [0,12800,65800,34100,46700,47500,82800,5480,735,578,0];

days_NK = [0,1,2,3,4,5,6,7,8,9,10, 11, 14, 17, 21, 28, 64, 67];

ki67_NK_ct = [23.2,15.8,15.1,34.8,34.7,48.6,67.5,122.3,23.6,79.1,42.3,...
              89.3,34.8,12.9,21.6,17.5,38.7,17.1];

tforward = (0:0.1:67);
tforward_v = (0:0.1:20);
tforward_n = (0:0.1:70);
tmeasure_v = [1:10:101];

tmeasure_n = [1 11 21 31 41 51 61 71 81 91 101 111 141 171 211 281 641 671];




function dy = model_Zika(t,y,k)

 dy = zeros(4,1); 

 % SM2
 dy(1) = -k(1)*y(3)*y(1); % S susceptible target cells
 dy(2) =  k(1)*y(3)*y(1) - k(2)*y(2)*y(4); % I infected target cells
 dy(3) =  y(2) - k(3)*y(3); % V free virus
 dy(4) =  k(4)*y(2)*y(2)*y(4)/(y(2)*y(2)+k(5)) - k(6)*y(4); % N natural killer cells

end



function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,[350000 1 0 23]);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  BN = y(tmeasure_n(:),4)';
  BV = y(tmeasure_v(:),3)'; 
  
 error_in_data = sum((BN - ki67_NK_ct).^2) + 0.00001*sum((BV - plasmaViral).^2);
 
 
end

%  k=1.0e+06 *[0.000000000021441   0.000000024945524   0.000002264335954   0.000004399153869...
%              1.444267002587110   0.000003002785340]; 

% k= 1.0e+05 *[0.000000000217478   0.000000400685651   0.000021460730788   0.000041739872158...
%                    6.594230516358982   0.000030590492118];
% k=[ 0.000247945206888   0.111794452799288   0.278702977167176   0.2849666079068...
%           0.000000000257592   0.061484469962943 ];
%      k=[  0.000641666882343   0.148070661543721   0.214754226699969   0.235447377298823...
%               0.000000000191736   0.028241921677374];
%         k=[  0.000644519418335   0.148178016062530   0.215781085566885   0.236814388571059...
%                   0.000000000170637   0.028085276618768 ];
%  

k = [0.000645284626980   0.148331268373208   0.215587884727577   0.236471262843213   0.000000000169701   0.027724940405331]            
 lb = [0 0 0 0 0 0];

 %[k,fval] = fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter', 'TolX', 1e-14,'TolFun',1e-14,...
 %                      'MaxIter',5000,'MaxFunEvals',5000)); 

 disp(k);


 [T,YV] = ode45(@(t,y)(model_Zika(t,y,k)),tforward_v,[350000 1 0 23]);
 [T,YN] = ode45(@(t,y)(model_Zika(t,y,k)),tforward_n,[350000 1 0 23]);
 
 
 
 figure(1)
 plot(days,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward_v, (YV(:,3)), '-b','LineWidth',3)

xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Male Monkey 393422'},'FontWeight','normal')
set(gca,'YTick',[0, 2*10^4, 4*10^4, 6*10^4, 8*10^4, 10*10^4],'YTickLabel',...
    {'0', '2x10^4', '4x10^4', '6x10^4', '8x10^4', '10^5'}, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')
 hold off
 
 
 figure(2)
 plot(days_NK,ki67_NK_ct,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward_n, (YN(:,4)), '-b','LineWidth',3)

ylabel('Amount of Activated Natural Killer Cells','FontSize',14,'FontName','Sans-serif');
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Activated Natural Killer cells in male monkey 393422'},'FontWeight','normal')
set(gca,'YLim',[0 150], 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
 hold off

end

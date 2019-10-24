function Fitting_Zika_Within_SM2_Monkey181856

clear all
close all
clc
format long

%MONKEY 181856 DATA FEMALE NONPREGNANT MONKEY

days = [0,1,2,3,4,5,6,7,8,9,10];

plasmaViral = [0,5800,168000,104000,26100,4220,2380,803,403,0,0];

days_NK = [0,2,4,6,7,8,9,10,14,21,28,64,70,71,74,77,80,84,91,99];

ki67_NK_ct = [13.3, 4.85, 11.79,26.34, 15.73, 40.17,59.14,21.85,21.31,...
              44.92,18.53,23.09,11.76,14.76,14.72,15.76,12.38,7.58,17.33, 17.34];

tforward = (0:0.1:100);
tforward_v = (0:0.1:10);
tforward_n = (0:0.1:99);
tmeasure_v = [1:10:101];

tmeasure_n = [1 21 41 61 71 81 91 101 141 211 281 641 701 711 741 771 801 841 911 991];

ic=[650000 1 0 13];


function dy = model_Zika(t,y,k)

 dy = zeros(4,1); 

 % SM2
 dy(1) = -k(1)*y(3)*y(1); % S susceptible target cells
 dy(2) =  k(1)*y(3)*y(1) - k(2)*y(2)*y(4); % I infected target cells
 dy(3) =  y(2) - k(3)*y(3); % V free virus
 dy(4) =  k(4)*y(2)*y(2)*y(4)/(y(2)*y(2)+k(5)) - k(6)*y(4); % N natural killer cells

end



function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,ic);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  BN = y(tmeasure_n(:),4)';
  BV = y(tmeasure_v(:),3)'; 
  
 error_in_data = sum((BN - ki67_NK_ct).^2) + 0.0000001*sum((BV - plasmaViral).^2);


end

%  k=1.0e+06 *[0.00000000021441   0.000000024945524   0.000002264335954   0.0000001399153869...
%              1.444267002587110   0.00000003002785340]; 

         
% k=1.0e+09 *[0.000000000000183   0.000000000063953   0.000000001719360   0.000000000526790...
%                   4.914364763258036   0.000000000010018];

k = [0.000221748110276   0.065794496289213   1.887523612110624   0.132039541037805...
    48.952546725733022   0.012624064846217]
%results in 
k =  [0.000221748108896   0.065794496805742   1.887523602280112...
      0.132039540684570  48.952530139771028   0.012624064866131];
  
  lb = [0 0 0 0 0 0];
% 
 %[k,fval] = fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter', 'TolX', 1e-14,'TolFun',1e-14,...
  %                     'MaxIter',5000,'MaxFunEvals',5000)); 

 disp(k);


 [T,YV] = ode45(@(t,y)(model_Zika(t,y,k)),tforward_v,ic);
 [T,YN] = ode45(@(t,y)(model_Zika(t,y,k)),tforward_n,ic);
 
 
 
 figure(1)
 plot(days,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward_v, (YV(:,3)), '-b','LineWidth',3)

xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus in non-pregnant female monkey 181856'},'FontWeight','normal')
ylim([0 2*10^5])
set(gca,'YTick',[0, 5*10^4, 10^5, 2*10^5],'YTickLabel',...
    {'0', '5x10^4', '10^5', '2x10^5'}, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')
 hold off
 
 
 figure(2)
 plot(days_NK,ki67_NK_ct,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
            hold on 
 plot(tforward_n, (YN(:,4)), '-b','LineWidth',3)

ylabel('Amount of Activated Natural Killer Cells','FontSize',14,'FontName','Sans-serif');
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Activated Natural Killer cells in non-pregnant female monkey 181856'},'FontWeight','normal')
set(gca,'YLim',[0 70], 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
 hold off

end

function k=Fitting_Zika_within_Scaled_SM3_1stTrimester

clear all
close all
clc
format long

%Pregnant Monkey 827577

days_v = [0,1,2,3,4,5,7,8,9,10,14,21,28,35,43,50,57,64,71,78];

plasmaViral = [0,2630,413000,668000,111000,7470,2210,8380,3770,...
               2550,512,3180,11500,3780,268,715,268,791,167,0];
days_n = [0,2,3,7,8,9,10,14,21,28,35,42,57,78];
ki67_NK_ct = [27,56,31,81,89,55,55,36,21,26,12,6,4,6];
     
tforward = (0:0.1:278);

tmeasure_v = [1,11,21,31,41,51,71,81,91,101,141,211,281,351,431,501,571,641,711,781];
tmeasure_n = [1,21,31,71,81,91,101,141,211,281,351,421,571,781];
     


ic = [45000000 10 0 28 1000000 0 0]; 




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
  BV = y(tmeasure_v(:),3)'; 
  
 error_in_data = sum((BN - ki67_NK_ct).^2) + 0.00001*sum((BV - plasmaViral).^2);
 

end

 
 % [beta   delta    pi     p_f, ... 
 % c       alpha    gamma  d  ...
 % beta_f  delta_f  p      pi_f ...
 % c_f]


%  k= 1.0e+02 *[ 0.00000124469023   0.000322268082967   0.000085404671   0.000000085404671 ...
%                0.023484331514236   0.000543965345304   7.762215836432640   0.000283197696364 ...
%                0.0100000002926   0.010000000000606   8.801040535883160   0.010001000000606 ...
%                0.000000000002926  ];
%  k= 1.0e+02 *[ 0.000000824469023   0.000322268082967   0.000085404671   0.000000085404671 ...
%                0.023484331514236   0.000543965345304   7.762215836432640   0.000283197696364 ...
%                0.0100000002926   0.010000000000606   8.801040535883160   0.010001000000606 ...
%                0.000000000002926  ]; %best start. Run 2000 iterations.
%  k= 1.0e+02 *[0.000000131225047   0.000544970851812   0.001720205731146   0.000000121075483   0.029532433747920 ...
%      0.001174971257014   0.022141685709209   0.000128984827969   0.000027318619907   0.019200820490006 ...
%            4.240122208525726   0.013092292012491   0.000000000033494];
%   k= 1.0e+02 *[ 0.000000128542956   0.000582390672182   0.001748564329308   0.000000034049151   0.026365498963020 ...
%       0.002104063669931   0.009692245498392   0.000772495369914   0.000005182771654   0.000485241133230 ...
%       4.266355893553571   0.035875671414473   0.000000000107090];
%   k= 1.0e+04 *[ 0.000000001293880   0.000005453672613   0.000017466605373   0.000000000042553   0.000274234785426 ...
%           0.000024059479408   0.001994637615029   0.000009473271378   0.000000124196689   0.000002901563981 ...
%           1.031806719391026   0.000373428082921   0.000000000002070];
%   k= 1.0e+03 * [0.000000012809299   0.000052622370104     0.002841106527057 ...
%                 0.000245626878444   0.564470154293374     0.000086475636612...
%                 0.000003360937923   0.000082760546604     0.00005112746877378540     0.000000000085409];
%             
% k = 1.0e+03 *[0.000000006596724   0.000028229718304   0.029314841452244...
%               0.000418710434924   6.623789097121285   0.000231498417149...
%               0.000000118658375   0.001375859591521   0.000105232339539...
%               0.000000001012118];
% 
% k = 1.0e+05 *[0.000000000050596   0.000000415076032   0.000199240144656...
%               0.000003673841816   3.496540954820121   0.000001690830073...
%               0.000000000000000   0.000292687319529   0.000000325503115...
%               0.000000000005258];
          
% k = 1.0e+06 *[0.000000000005042   0.000000042747356   0.000019812655309...
%               0.000000480624163   0.089067923037223   0.000000282560795...
%               0.000000000012200   3.575532474414540   0.000000044219153...
%               0.000000021203511];
% k = 1.0e+06 *[0.000000000005048   0.000000042511886   0.000019891257371...
%               0.000000484680417   0.144985290777617   0.000000278468889...
%               0.000000000006484   3.502939469909511   0.000000054918535...
%               0.000000021777292]; 
          
% k = 1.0e+06 *[0.000000000005050   0.000000042582065   0.000019896578914...
%               0.000000485150920   0.147493417946678   0.000000279266500...
%               0.000000000006985   3.371075105374305   0.000000056194507...
%               0.000000022036893]; 
          
k = 1.0e+06 *[0.000000000005050   0.000000042603184   0.000019898046524...
              0.000000485060891   0.147716353104356   0.000000279383936...
              0.000000000007006   3.504265363595249   0.000000056280792...
              0.000000022048771];
          
lb = [0 0 0 0 0 0 1e-8 0 0 0];
                  

%[k,~] = fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxIter',5000,'MaxFunEvals',5000))

 [T,Y] = ode23s(@(t,y)(model_Zika(t,y,k)),tforward,ic);
 
 
 
 %figure1 = figure('position', [0, 0, 1200, 400]); 
 figure(1)
 plot(days_v,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
 plot(tforward, (Y(:,3)), '-b','LineWidth',3)
axis([0 80 0 10^7])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Pregnant Monkey 827577','ZIKV Infected At The 1st Trimester'},'FontWeight','normal')
set(gca,'yscale','log','YTick',[10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7],'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')
 hold off
 
 %figure2 = figure('position', [0, 0, 1200, 400]); 
 figure(2)
 plot(days_n,ki67_NK_ct,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
  plot(tforward, (Y(:,4)), '-b','LineWidth',3)
axis([0 80 0 100])
ylabel('Natural Killer Cells','FontSize',14,'FontName','Sans-serif');
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
set(gca, 'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
title ({'Natural Killer Cell Counts In Pregnant Monkey 827577'},'FontWeight','normal')

 hold off
 
 figure(3)
 plot(days_v,plasmaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
 plot(tforward, (Y(:,3)), '-b','LineWidth',3)
%axis([0 80 0 10^4])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Pregnant Monkey 827577','ZIKV Infected At The 1st Trimester'},'FontWeight','normal')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')

 figure(4)

plot(tforward, (Y(:,7)), '-b','LineWidth',3)
%axis([0 80 0 10^4])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Fetus','Pregnant Monkey is Infected At The 1st Trimester'},'FontWeight','normal')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')


 figure(5)

plot(tforward(1:151), (Y(1:151,7)), '-b','LineWidth',3)
%axis([0 80 0 10^4])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Fetus','Pregnant Monkey is Infected At The 1st Trimester'},'FontWeight','normal')
set(gca,'yscale','log','YTick',[0 10^1, 10^2, 10^3, 10^4, 10^5],'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')

figure(6)
plot(tforward, (Y(:,3)), '-b','LineWidth',3)
hold on 
plot(tforward, (Y(:,7)), '-r','LineWidth',3)

xlabel('Days','FontSize',14,'FontName','Sans-serif' );
title ({'Free Virus In Pregnant Monkey 827577','ZIKV Infected At The 1st Trimester'},'FontWeight','normal')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Viral RNA copies per milliliter of plasma','FontSize',14,'FontName','Sans-serif')
end
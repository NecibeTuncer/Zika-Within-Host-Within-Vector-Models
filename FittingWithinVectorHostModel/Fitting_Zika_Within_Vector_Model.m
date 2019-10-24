function k= Fitting_Zika_Within_Vector_Model

clear all
close all
clc
format long


days_v = [2,3,4,5,6,7,10,14];

midgutViral = [3.89665, 3.97, 5.14, 5.574977, 5.38589, 5.5399, 5.756262, 5.33324];
days_n = [5, 6, 7, 10, 14];

salivaViral = [2.487, 4.0495, 3.5727, 6.4888, 6.785];
     
tforward = (0:0.1:14);

tmeasure_v = [21,31,41,51,61,71,101,141];
tmeasure_n = [51,61,71,101,141];
     


ic = [1000 0.0]; 




function dy = model_Zika(t,y,k)

 dy = zeros(2,1); 

 % Model 1
 dy(1) = k(1)*y(1)*(1-y(1)/k(2))-k(3)*y(1); % virus in the midgut
 dy(2) = k(3)*y(1)*y(1)/(k(4) + y(1)*y(1))+ k(5)*y(2)*(1-y(2)/k(6)); %virus in the salivary gland

 
end


function error_in_data = err_in_data(k)


 [t,y] = ode23s(@(t,y)(model_Zika(t,y,k)),tforward,ic);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  VM = log10(y(tmeasure_v(:),1)');
 % VH = y(tmeasure_v(:),2)'; 
  VS = log10(y(tmeasure_n(:),2)');
  
 error_in_data = sum((VM - midgutViral).^2) +  sum((VS - salivaViral).^2);
 

end


  

%  k =    1.0e+06 *[0.000001150607772   0.395179196949023   0.000000000076623...
%                   0.000001425417144    0.000001763870021    6.846296490949689];
              
 %results in 
%  k = 1.0e+06 *[0.000001258872224   0.432705216819770   0.000000108754431...   
%                0.002954458426231   0.000001769318844   6.838637317291384];
%            
%  
%  k = 1.0e+06 *[0.000001260828899   0.433377378612009   0.000000110712308...
%                0.019541590122065   0.000001767541515   6.843581891443305];
%            
%  k = 1.0e+06 *[0.000001260830998   0.433380458505429   0.000000110714319...
%                0.019533227991548   0.000001767537909   6.843605178932394];
           
%  k = 1.0e+06 *[0.000001260830433   0.433380850773171   0.000000110713399...
%                0.019529867237489   0.000001767538722   6.843610335784404];
 %results in           
 k = [1.26082924880034,433378.838892983,0.110713342233325,...
     19534.2728313170,1.76753893630032,6843584.00734574];          
 lb = [0 0 0  0.01 0 0];
                  
 % %for i=1:5
 %     i
 %[k,~] = fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000));

 disp(k);
%  end
 [T,Y] = ode23s(@(t,y)(model_Zika(t,y,k)),tforward,ic);
 
 
 
 figure(1)% = figure('position', [0, 0, 1200, 1400]);
 
 plot(days_v,midgutViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
 plot(tforward, (log10(Y(:,1))), '-b','LineWidth',3)
%axis([0 80 0 10^6])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
%title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
%set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
 %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Virus titer Log10 TCID50/ml in the midgut','FontSize',14,'FontName','Sans-serif')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
 hold off
 
 
 figure(2)% = figure('position', [0, 0, 1200, 1400])
 
 plot(days_n,salivaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
 plot(tforward, (log10(Y(:,2))), '-b','LineWidth',3)
  ylabel('Virus titer Log10 TCID50/ml in the saliva glands','FontSize',14,'FontName','Sans-serif');
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
 hold off 
 
 figure(3)
 plot(days_v,10.^midgutViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on
 plot(tforward, Y(:,1), '-b','LineWidth',3)
%axis([0 80 0 10^6])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
%title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
%set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
 %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Virus titer  TCID50/ml in the midgut','FontSize',14,'FontName','Sans-serif')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
  
 figure(4)
  plot(days_n,10.^salivaViral,'Marker','.','Color',[1 0 0],...
                'MarkerSize',30,'LineStyle','none')
             hold on 
 plot(tforward, Y(:,2), '-b','LineWidth',3)
%axis([0 80 0 10^6])
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
%title ({'Free Virus in the pregnant monkey 827577','ZIKV infected at the 1st Trimester'},'FontWeight','normal')
%set(gca,'YTick',[0,10^0,10^1,10^2,10^3,10^4,10^5,10^6],'YTickLabel',...
 %  {'0','10^0','10^1','10^2','10^3','10^4','10^5','10^6'},'yscale','log','LineWidth',2,'FontSize',14,'FontName','Sans-serif');
ylabel('Virus titer  TCID50/ml in the saliva glands','FontSize',14,'FontName','Sans-serif')
set(gca,'LineWidth',2,'FontSize',14,'FontName','Sans-serif');
  
  
end
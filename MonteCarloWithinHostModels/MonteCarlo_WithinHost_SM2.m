function  total_ARE = MonteCarlo_WithinHost_SM2

%Monkey181856

clear all
close all
clc

numiter = 1000; 


X = zeros(6,numiter);


tforward = (0:0.1:100);

tmeasure_v = [1:10:101];

tmeasure_n = [1 21 41 61 71 81 91 101 141 211 281 641 701 711 741 771 801 841 911 991];

ic=[650000 1 0 13];

true_params =  [0.000221748108896   0.065794496805742   1.887523602280112...
                0.132039540684570  48.952530139771028   0.012624064866131];


[~,y_trp] = ode45(@(t,y)(model_Zika(t,y,true_params)),tforward,ic);



function dy = model_Zika(t,y,k)

 dy = zeros(4,1); 

 % SM2
 dy(1) = -k(1)*y(3)*y(1); %S susceptible target cells
 dy(2) =  k(1)*y(3)*y(1) - k(2)*y(2)*y(4); %I infected target cells
 dy(3) =  y(2) - k(3)*y(3); %V free virus
 dy(4) =  k(4)*y(2)*y(2)*y(4)/(y(2)*y(2)+k(5)) - k(6)*y(4); %N natural killer cells

end



function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)(model_Zika(t,y,k)),tforward,ic);%,...;
                  %odeset('RelTol',10^(-2),'AbsTol',10^(-2),'NormControl','off'));

  BN = y(tmeasure_n(:),4)';
  ZV = y(tmeasure_v(:),3)'; 
  
 error_in_data = sum((BN - ydata_n).^2) + 0.0000001*sum((ZV - ydata_v).^2);


end

 

noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:length(noiselevel)
  
rng default
noiselev = noiselevel(noisei)

for i= 1:numiter
    i

   ytrp_v = y_trp(tmeasure_v(:),3);
   ydata_v = noiselev*randn(1,size(tmeasure_v,2)) + ytrp_v';
   
   ytrp_n = y_trp(tmeasure_n(:),4);
   ydata_n = noiselev*randn(1,size(tmeasure_n,2)) + ytrp_n';

k = true_params;

 
 lb = [0 0 0 0 0 0];
 
k =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter')); 
                        
X(:,i) = k';


end


arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;
 
end

save('MCS_SM2_ARE')

end
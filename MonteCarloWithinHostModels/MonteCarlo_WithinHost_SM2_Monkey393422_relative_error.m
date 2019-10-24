function  total_ARE = MonteCarlo_WithinHost_SM2_Monkey393422_relative_error

%Monkey393422

clear all
close all
clc

numiter = 1000; 


X = zeros(6,numiter);


tforward = (0:0.1:67);

tmeasure_v = [1:10:101];

tmeasure_n = [1 11 21 31 41 51 61 71 81 91 101 111 141 171 211 281 641 671];

ic=[350000 1 0 23];

true_params = [0.000645284626980   0.148331268373208   0.215587884727577...
               0.236471262843213   0.000000000169701   0.027724940405331];   


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
  
 error_in_data = sum((BN - ydata_n).^2) + 0.00001*sum((ZV - ydata_v).^2);


end

 

noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:length(noiselevel)
  
rng default
noiselev = noiselevel(noisei)

for i= 1:numiter
    i

   ydata_v = normrnd(y_trp(tmeasure_v(:),3)',noiselev*y_trp(tmeasure_v(:),3)');
   
   ydata_n = normrnd(y_trp(tmeasure_n(:),4)',noiselev*y_trp(tmeasure_n(:),4)');

k = true_params;

 
 lb = [0 0 0 0 0 0];
 
k =  fminsearchbnd(@err_in_data,k,lb,[]); 
                        
X(:,i) = k';


end


arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;
 
end
save('MCS_SM2_Monkey_393422_ARE_relative_error')
end
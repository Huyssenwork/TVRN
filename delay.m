%import data
Pomega=xlsread('cov_19_gen.xlsx','B2:B10');
data=xlsread('New Zealand.xlsx','B2:B301');
data_size=size(data);
I_pre=[data(1)];
%serial interval distribution 
nday=data_size(1);
nul=zeros(nday, 1);
Pomega=cat(1,Pomega,nul);
Rgrid = linspace(0.1, 3,1000);
%delay distribution 
Gp=[0.15,0.26,0.23,0.16,0.12,0.05,0.03]';
Gp=cat(1,Gp,nul);
% Prob vector for R and prior
pR = zeros(nday, 1000); 
pR(1, :) = 0.001;
% Mean, median and confidence on R
Rm = zeros(1, nday); 
% Initial stats
Rm(1) = pR(1, :)*Rgrid'; 
% Precompute state distributions
pstate = zeros(1000, 1000);
for j = 1:1000
    pstate(j, :) = normpdf(Rgrid(j), Rgrid, sqrt(Rgrid)*0.01);
end
Lam=zeros(nday,1);
% Update prior to posterior sequentially
Dday_es=[10];
I_pre=[data(1)];
for i = 2:nday
    %计算延迟
    Gpat = Gp(2:i+5);
    m = sum(Dday_es(i+6:-1:3)'.*Gpat); 
    %Dday_es=[Dday_es data(i)-m];
    %计算总传染性
    Pomegat = Pomega(1:i-1);
    Lam(i) = sum(data(i-1:-1:1).*Pomegat);
    % Compute rate from Poisson renewal
    rate = Lam(i).*Rgrid;
    
    rate =Gp(1).*rate;
    mu=rate+m;
    % Probabilities of observations
    pI = poisspdf(data(i), mu);
    % State equations for R
    pRup = pR(i-1, :)*pstate;
    % Posterior over R (updated)
    pR(i, :) = pRup.*pI;
 
    pR(i, :) = pR(i, :)/sum(pR(i, :));
    % Posterior mean and CDF
    Rm(i) = pR(i, :)*Rgrid';
    %Rm(i)=Rm(i)^(1/3);
    d=Rm(i)*Lam(i);
    
    Dday_es=[Dday_es d];
    I_pre=[I_pre,Gp(1)*d+m];
end



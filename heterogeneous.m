%import data
Pomega=xlsread('cov_19_gen.xlsx','B2:B10');
data=xlsread('New Zealand.xlsx','B2:B301');
R_len=200;
Rgrid = linspace(0.1,3,R_len);
r_len=10;
rgrid=linspace(0.4,0.6,r_len);
% Prob vector for R and r and prior
pR_r=zeros(nday*R_len,r_len); 
pR_r_day=R_len*ones(1,nday);
pR_r=mat2cell(pR_r,pR_r_day,[r_len]);
a=1/(R_len*r_len);
pR_r{1,1} = (1/(R_len*r_len))*ones(R_len,r_len);
% Mean, median and confidence on R
Rm = zeros(1, nday); 
rm = zeros(1, nday); 
% Initial stats
pR=sum(pR_r{1,1},2);
Rm(1) = pR'*Rgrid'; 
pr=sum(pR_r{1,1},1);
rm(1) = pr*rgrid';
% stata distribution 
pstate=ones(R_len*R_len,r_len*r_len);
lenth=R_len*ones(1,R_len);
lenth2=r_len*ones(1,r_len);
pstate=mat2cell(pstate,lenth,lenth2);
for m=1:R_len
    for n=1:r_len
        cell_matrix=ones(R_len,r_len);
        for j=1:R_len
            for i=1:r_len
                cell_matrix(j,i)=mvnpdf([Rgrid(m),rgrid(n)],[Rgrid(j),rgrid(i)],[0.05*Rgrid(j),0;0,0.1*rgrid(i)]);
            end
        end
        pstate{m,n}=cell_matrix;
    end
end

% Update prior to posterior sequentially
for i = 2:nday  
    Pomegat = Pomega(1:i-1);
    Lam(i) = sum(data(i-1:-1:1).*Pomegat);
    % state equation
    pRup=pstate;
    for k=1:R_len
        for l=1:r_len
           pRup{k,l}=pstate{k,l}.*pR_r{i-1,1};
        end
    end
    pRup=cellfun(@(x) sum(x(:)),pRup);
    % Compute rate from Poisson renewal
    rate = Lam(i).*Rgrid;
    %rate =Gp(1).*rate;
    mu=rate;
    pI=ones(R_len,r_len);
    for j=1:r_len
         mu=mu+rgrid(j);
         p=rgrid(j).*(mu.^-1);
        pI(:,j) = nbinpdf(data(i),rgrid(j),p);
    end
    % Posterior over R (updated)
    pR_r{i,1} = pRup.*pI;
    pR=sum(pR_r{i,1},2);
    pR=pR/sum(pR);
    pr=sum(pR_r{i,1},1);
    pr=pr/sum(pr');
    pR_r{i,1}=pR*pr;
    % Posterior mean and CDF
    Rm(i) = pR'*Rgrid';
    rm(i) = pr*rgrid';
    d=Rm(i)*Lam(i);
    
end

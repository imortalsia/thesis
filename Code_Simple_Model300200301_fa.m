%%                 Section 1: Plot Wealth Inequality
clear; clc;
[data1,text1,~]  = xlsread('country_data','WID_top','A1:E100');
[data2,text2,~]  = xlsread('country_data','WID_gini','A1:E100');
country_Top_WID  = text1(1,2:end);
country_Gini_WID = text2(1,2:end);
Year_Top_WID  = data1(:,1);
Year_Gini_WID = data2(:,1);
Top_WID  = data1(:,2:end);
Gini_WID = data2(:,2:end);
color_Top_WID  = [1 0 0; .39216 .58431 .92941; .25098 .87843 .81569; 0 0 0];
color_Gini_WID = [1 0 0; .39216 .58431 .92941; .25098 .87843 .81569; 0 0 0];
figure(1);
subplot(1,2,1);
    for i=1:size(country_Top_WID,2)
        plot(Year_Top_WID,Top_WID(:,i),'Color',color_Top_WID(i,:),'LineWidth',2,'Marker','.','MarkerSize',20);
        hold on;
    end
    hold off;
    axis([Year_Top_WID(1) Year_Top_WID(end) 0.4 0.8]);
    set(gca,'ytick',0.4:0.1:0.8);
    set(gca,'yticklabel',{'40%','50%','60%','70%','80%'},'FontSize',12,'FontName','Times New Roman');
    title('(a) Top10% Wealth Share','FontSize',12,'FontName','Times New Roman')
    L1=legend(country_Top_WID,'location','SouthOutside','Orientation','horizontal','Box','off','FontSize',12,'FontName','Times New Roman');
    L1.Position = [0.455 0.03 0.1 0.03];
subplot(1,2,2);
    for i=1:size(country_Gini_WID,2)
        plot(Year_Gini_WID,Gini_WID(:,i),'Color',color_Gini_WID(i,:),'LineWidth',2,'Marker','.','MarkerSize',20);
        hold on;
    end
    hold off;
    axis([Year_Gini_WID(1) Year_Gini_WID(end) 0.5 0.9]);
    set(gca,'ytick',0.5:0.1:0.9,'FontSize',12);
    title('(b) Wealth Gini coefficient','FontSize',12,'FontName','Times New Roman')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 10 6]);
print(gcf,'-dpng','-r1080','Figure1.Wealth Inequality in Selected Countries.jpg');

%%                 Section 2: Parameter setting
clear; clc;
npel = 300;          % Number of agents
nper = 300;          % Number of periods
W0 = ones(npel,1)*3;  % Initial wealth endowment
Y  = npel*1;          % Aggregat output in each period
mid_per= nper*0.3;    % The starting period for social mobility 30 years
top    = 0.1;         % Top 10% share
bottom = 0.5;         % Bottom 50% share
TM_top = 0.1;         % Top percentage for the transition matrix
TM_btm = 0.5;         % Bottom percentage for the transition matrix
C0     = 0.4;         % Subsistence consumption over GDP per capita
C1     = 0.6;         % Consumption propensity
beta1  = 1;           % Human capital return
sigma  = 0.1;         % Dispersion of (log) human capital return
growth = 0.08;        % Net annual growth rate
tax    = 0.2;         % Proportional income tax rate
allowance = 1;        % Tax allowance relative to income
upper  = 0.5;         % +50% variation of the key parameter
lower  = 0.5;         % -50% variation of the key parameter
G_upper= 0.16;        % Growth rate variated by +100%
G_lower= 0.008;       % Growth rate variated by -90%
Mob_dim= 5;           % Dimension of the social mobility matrix (= 5 or 10)
Lorenz_per   = nper;  % set the period at which we plot the Lorenz curve
W_sorted_end0= zeros(npel,4);
W_sorted_end1= zeros(npel,4);
mob_mat = zeros(Mob_dim,Mob_dim,5);   % Social Mobility matrix (mid to end)
Shorrocks_S  = zeros(5,1);            % Static Shorrocks Index
save settings_simple.mat;             % All settings for later import

%%                 Section 3: Introductory experiment
clear; clc;
load settings_simple.mat;
[sim_wealth,~] = introduction(W0,Y,nper,C0,C1);
sort_wealth = sort(sim_wealth);
ylow = floor(min(sim_wealth(:,nper)));
yup  = fix(max(sim_wealth(:,nper)));
yaxis_low = ylow*1.1;
yaxis_up  = yup*1.1;
npel_aix  = 1:npel;
figure(2);
subplot(1,2,2);
   set(gca,'Position',[0.5,0.15,0.4,0.75]);
   plot(npel_aix,W0','k--','LineWidth',1);                         hold on;
   plot(npel_aix,sort_wealth(:,nper),'k','LineWidth',1);           hold on;
   text(533,18,'Final','FontName','Times');
   arrow1 = annotation('arrow','headStyle','cback1','HeadLength',5,'HeadWidth',5);
   set(arrow1,'position',[0.74 0.60 0 -0.05]);                     hold on;
   text(533,-8,'Initial','FontName','Times');
   arrow2 = annotation('arrow','headStyle','cback1','HeadLength',5,'HeadWidth',5);
   set(arrow2,'position',[0.74 0.43 0 0.06]);                      hold off;
   axis_value = axis;
   axis([axis_value(1) axis_value(2) yaxis_low yaxis_up]);
   set(gca,'ytick',ylow:(mean(W0)-ylow)/2:yup);
   set(gca,'YAxisLocation','right');
   set(gca,'yTicklabel',{'','','','','','',});
   w_position = yup*1.3;
   ylabel('Wealth','Position',[0,w_position],'Rotation',0,'FontName','Times New Roman');
   xlabel('Individual (sorted)','FontName','Times New Roman');
subplot(1,2,1); % Don't use 'histogram' as hist is inconsistent with the pdf
   set(gca,'Position',[0.1,0.15,0.4,0.75]);
   [fp,xp]=ecdf(sim_wealth(:,nper));
   ppx=min(xp):(max(xp)-min(xp))/(50-1):max(xp);
   ppy=ecdfhist(fp,xp,50);
   ppy2=ppy;
   barh(ppx,ppy2,'w');                                             hold on;
   x_pdf=sort_wealth(:,nper);
   y_pdf=normpdf(x_pdf,mean(W0),sqrt(nper*(npel-1)/npel)*(1-C1));
   plot(y_pdf,x_pdf,'k','LineWidth',1);                            hold on;
   text(0.035,30,'Asymptotically Normal','FontName','Times');
   arrow3 = annotation('arrow','headStyle','cback1','HeadLength',5,'HeadWidth',5);
   set(arrow3,'position',[0.22 0.7 0.06 -0.06]);                   hold off;
   clear axis_value;
   axis_value = axis;    
   axis([axis_value(1) axis_value(2) yaxis_low yaxis_up]);
   set(gca,'XAxisLocation','bottom');
   set(gca,'XDir','reverse');   % turn A axis
   set(gca,'XTicklabel',{'',0.01,0.02,0.03,0.04});
   set(gca,'ytick',ylow:(mean(W0)-ylow)/2:yup);
   set(gca,'YAxisLocation','right');
   xlabel('Probability','FontName','Times New Roman');
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 6 4.5]);
print(gcf,'-dpng','-r1080','Figure2.Introductory experiment.jpg');

%%                 Section 4.1: Baseline Model
clear; clc;
load settings_simple.mat;
[sim_W,~,~] = fit_sim(W0,0,0,0,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
sort_W = sort(sim_W);
%------Final top and bottom shares
top_share_base=sum(sort_W(npel*(1-top)+1:npel,nper),1)/sum(sort_W(:,nper),1);
btm_share_base=sum(sort_W(1:npel*bottom,nper),1)/sum(sort_W(:,nper),1);
%------Lorenz Curve & Gini (Don't scale)
acc_W_base  = cumsum(sort_W(:,Lorenz_per))/sum(sort_W(:,Lorenz_per));
W_Gini_base = gini_coeff(acc_W_base);
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_W(npel*(1-TM_top),nper); % Lower bound of the top "TM_top"
Threshold_btm = sort_W(npel*TM_btm,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
up_rate_base  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate_base= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<=Threshold_btm)/(npel*TM_top);
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=1;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
Shorrocks_S(k,1) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index 
%--------------------------------------------------------------------------
%                 Section 4.2: Income Tax Model
tax_upper = tax*(1+upper);
tax_lower = tax*(1-lower);
[sim_W,~,~] = fit_sim(W0,0,0,1,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
[sim_W_upper,~,~] = fit_sim(W0,0,0,1,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax_upper,allowance);
[sim_W_lower,~,~] = fit_sim(W0,0,0,1,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax_lower,allowance);
sort_W = sort(sim_W);
sort_W_upper = sort(sim_W_upper);
sort_W_lower = sort(sim_W_lower);
% Final top and bottom shares
top_share_tax=sum(sort_W(npel*(1-top)+1:npel,nper),1)/sum(sort_W(:,nper),1);
btm_share_tax=sum(sort_W(1:npel*bottom,nper),1)/sum(sort_W(:,nper),1);
%------Lorenz Curve & Gini (Don't scale)
acc_W_tax = cumsum(sort_W(:,Lorenz_per))/sum(sort_W(:,Lorenz_per));
acc_W_tax_upper = cumsum(sort_W_upper(:,Lorenz_per))/sum(sort_W_upper(:,Lorenz_per));
acc_W_tax_lower = cumsum(sort_W_lower(:,Lorenz_per))/sum(sort_W_lower(:,Lorenz_per));
W_Gini_tax = gini_coeff(acc_W_tax);
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_W(npel*(1-TM_top),nper); % Lower bound of the top "TM_top"
Threshold_btm = sort_W(npel*TM_btm,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
up_rate_tax  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate_tax= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<=Threshold_btm)/(npel*TM_top);
W_sorted_end0(:,1) = sort_W(:,nper);         % Save the end-period wealth
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=2;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
Shorrocks_S(k,1) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index
%--------------------------------------------------------------------------
%                 Section 4.3: Constant Growth Model
[sim_W,~,~] = fit_sim(W0,0,1,0,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
[sim_W_upper,~,~] = fit_sim(W0,0,1,0,0,0,Y,nper,C0,C1,beta1,sigma,G_upper,tax,allowance);
[sim_W_lower,~,~] = fit_sim(W0,0,1,0,0,0,Y,nper,C0,C1,beta1,sigma,G_lower,tax,allowance);
sort_W = sort(sim_W);
sort_W_upper = sort(sim_W_upper);
sort_W_lower = sort(sim_W_lower);
% Final top and bottom shares
top_share_growth = sum(sort_W(npel*(1-top)+1:npel,nper),1)/sum(sort_W(:,nper),1);
btm_share_growth = sum(sort_W(1:npel*bottom,nper),1)/sum(sort_W(:,nper),1);
%------Lorenz Curve & Gini (Don't scale)
acc_W_growth = cumsum(sort_W(:,Lorenz_per))/sum(sort_W(:,Lorenz_per));
acc_W_growth_upper = cumsum(sort_W_upper(:,Lorenz_per))/sum(sort_W_upper(:,Lorenz_per));
acc_W_growth_lower = cumsum(sort_W_lower(:,Lorenz_per))/sum(sort_W_lower(:,Lorenz_per));
W_Gini_growth = gini_coeff(acc_W_growth);
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_W(npel*(1-TM_top),nper); % Lower bound of the top "TM_top"
Threshold_btm = sort_W(npel*TM_btm,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
up_rate_growth  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate_growth= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<=Threshold_btm)/(npel*TM_top);
W_sorted_end0(:,2) = sort_W(:,nper);
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=3;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
Shorrocks_S(k,1) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index 
%--------------------------------------------------------------------------
%                 Section 4.4: Human Capital Model
beta1_upper = beta1*(1+upper);
beta1_lower = beta1*(1-lower);
[sim_W,~,~]  = fit_sim(W0,1,0,0,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
[sim_W_upper,~,~] = fit_sim(W0,1,0,0,0,0,Y,nper,C0,C1,beta1_upper,sigma,growth,tax,allowance);
[sim_W_lower,~,~] = fit_sim(W0,1,0,0,0,0,Y,nper,C0,C1,beta1_lower,sigma,growth,tax,allowance);
sort_W = sort(sim_W);
sort_W_upper = sort(sim_W_upper);
sort_W_lower = sort(sim_W_lower);
% Final top and bottom shares
top_share_effort = sum(sort_W(npel*(1-top)+1:npel,nper),1)/sum(sort_W(:,nper),1);
btm_share_effort = sum(sort_W(1:npel*bottom,nper),1)/sum(sort_W(:,nper),1);
%------Lorenz Curve & Gini (Don't scale)
acc_W_human = cumsum(sort_W(:,Lorenz_per))/sum(sort_W(:,Lorenz_per));
acc_W_human_upper = cumsum(sort_W_upper(:,Lorenz_per))/sum(sort_W_upper(:,Lorenz_per));
acc_W_human_lower = cumsum(sort_W_lower(:,Lorenz_per))/sum(sort_W_lower(:,Lorenz_per));
W_Gini_human = gini_coeff(acc_W_human);
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_W(npel*(1-TM_top),nper); % Lower bound of the top "TM_top"
Threshold_btm = sort_W(npel*TM_btm,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
up_rate_human  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate_human= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<=Threshold_btm)/(npel*TM_top);
W_sorted_end0(:,3) = sort_W(:,nper);
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=4;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
Shorrocks_S(k,1) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index 
%--------------------------------------------------------------------------
%                 Section 4.5: Physical Capital Model
% As no investment rate is involved, do not consider the upper & lower bounds
[sim_W,~,~] = fit_sim(W0,0,0,0,1,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
sort_W = sort(sim_W);
% Final top and bottom shares
top_share_invest = sum(sort_W(npel*(1-top)+1:npel,nper),1)/sum(sort_W(:,nper),1);
btm_share_invest = sum(sort_W(1:npel*bottom,nper),1)/sum(sort_W(:,nper),1);
%------Lorenz Curve & Gini (Don't scale)
acc_W_invest = cumsum(sort_W(:,Lorenz_per))/sum(sort_W(:,Lorenz_per));
W_Gini_capital = gini_coeff(acc_W_invest);
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_W(npel*(1-TM_top),nper); % Lower bound of the top "TM_top"
Threshold_btm = sort_W(npel*TM_btm,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
up_rate_invest  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate_invest= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<=Threshold_btm)/(npel*TM_top);
W_sorted_end0(:,4) = sort_W(:,nper);
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=5;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
Shorrocks_S(k,1) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index

%%                 Section 5: Plot Thought Experiments
xForFill   = [linspace(0,1,npel),fliplr(linspace(0,1,npel))];
yForFill_1 = [acc_W_tax_lower',fliplr(acc_W_tax_upper')];
yForFill_2 = [acc_W_growth_lower',fliplr(acc_W_growth_upper')];
yForFill_3 = [acc_W_human_lower',fliplr(acc_W_human_upper')];
figure(3)
subplot(2,2,1);
   set(gca,'Position',[0.1,0.57,0.37,0.37]);
   plot(linspace(0,1,npel),acc_W_base,'k--','LineWidth',0.5);      hold on;
   plot(linspace(0,1,npel),linspace(0,1,npel),'k','LineWidth',1);  hold on;
   plot(linspace(0,1,npel),acc_W_tax,'k','LineWidth',1.0);         hold on;
   fill(xForFill,yForFill_1,[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',[0.5 0.5 0.5]); % fill gray 
   hold off;
   xlabel('(a) Proportional income tax','FontSize',12,'FontName','Times New Roman');
subplot(2,2,2);
   set(gca,'Position',[0.57,0.57,0.37,0.37]);
   plot(linspace(0,1,npel),acc_W_base,'k--','LineWidth',0.5);      hold on;
   plot(linspace(0,1,npel),linspace(0,1,npel),'k','LineWidth',1);  hold on;
   plot(linspace(0,1,npel),acc_W_growth,'k','LineWidth',1.0);      hold on;
   fill(xForFill,yForFill_2,[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',[0.5 0.5 0.5]);
   hold off;
   xlabel('(b) Economic growth','FontSize',12,'FontName','Times New Roman');
subplot(2,2,3);
   set(gca,'Position',[0.1,0.1,0.37,0.37]);
   plot(linspace(0,1,npel),acc_W_base,'k--','LineWidth',0.5);      hold on;
   plot(linspace(0,1,npel),linspace(0,1,npel),'k','LineWidth',1);  hold on;
   plot(linspace(0,1,npel),acc_W_human,'k','LineWidth',1.0);       hold on;
   fill(xForFill,yForFill_3,[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',[0.5 0.5 0.5]);
   hold off;
   xlabel('(c) Human capital','FontSize',12,'FontName','Times New Roman');
subplot(2,2,4);
   set(gca,'Position',[0.57,0.1,0.37,0.37]);
   plot(linspace(0,1,npel),acc_W_base,'k--','LineWidth',0.5);      hold on;
   plot(linspace(0,1,npel),linspace(0,1,npel),'k','LineWidth',1);  hold on;
   plot(linspace(0,1,npel),acc_W_invest,'k','LineWidth',1.0);     hold off;
   xlabel('(d) Physical capital','FontSize',12,'FontName','Times New Roman');
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8 6]);
print(gcf,'-dpng','-r1080','Figure3.Lorenz Curves of Thought Experiments.jpg');

%%                 Section 6: Make Tables
% Table1: Wealth Inequality & Social Mobility of simple experiments
Wealth_Gini = [W_Gini_base;W_Gini_tax;W_Gini_growth;W_Gini_human;W_Gini_capital];
tops_all = [top_share_base;top_share_tax;top_share_growth;top_share_effort;top_share_invest];
btms_all = [btm_share_base;btm_share_tax;btm_share_growth;btm_share_effort;btm_share_invest];
up_rates   = [up_rate_base;up_rate_tax;up_rate_growth;up_rate_human;up_rate_invest];
down_rates = [down_rate_base;down_rate_tax;down_rate_growth;down_rate_human;down_rate_invest];
table1_value= num2cell([Wealth_Gini,tops_all,up_rates,down_rates,Shorrocks_S]);
table1_mix  = [{'Baseline';'Tax';'Growth';'Human Capital';'Physical Capital'},table1_value];
table1  = [{'Thought Experiment','Gini','Top 10%','Upward','Downward','SI'};table1_mix];
xlswrite('Generated Tables.xlsx',table1,'Table1.Inequality and Mobility');

% Table2a: Transition Matrix of the Baseline Model
table2_title= {'Initial(Row), Final(Column)','Bottom 20%','Lower 20%','Middle 20%','Upper 20%','Top 20%'};
table2  = [table2_title;{'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'},num2cell(mob_mat(:,:,1))];   
xlswrite('Generated Tables.xlsx',table2,'Table2.Baseline Transition');
%%
% Table10: Transition Matrices of simple Models
%---For Tax Model
table10a_mix  = [{'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'},num2cell(mob_mat(:,:,2))];
table10a  = [{'Tax Model','Bottom 20%','Lower 20%','Middle 20%','Upper 20%','Top 20%'};table10a_mix];   
xlswrite('Generated Tables.xlsx',table10a,'Table10.TE Transition Matrices','A1:F6');
%---For Growth Model
table10b_mix  = [{'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'},num2cell(mob_mat(:,:,3))];
table10b  = [{'Growth Model','Bottom 20%','Lower 20%','Middle 20%','Upper 20%','Top 20%'};table10b_mix];   
xlswrite('Generated Tables.xlsx',table10b,'Table10.TE Transition Matrices','A7:F12');
%---For Human Capital Model
table10c_mix  = [{'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'},num2cell(mob_mat(:,:,4))];
table10c  = [{'Human Capital Model','Bottom 20%','Lower 20%','Middle 20%','Upper 20%','Top 20%'};table10c_mix];   
xlswrite('Generated Tables.xlsx',table10c,'Table10.TE Transition Matrices','A13:F18');
%---For Physical Capital Model
table10d_mix  = [{'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'},num2cell(mob_mat(:,:,5))];
table10d  = [{'Physical Capital Model','Bottom 20%','Lower 20%','Middle 20%','Upper 20%','Top 20%'};table10d_mix];   
xlswrite('Generated Tables.xlsx',table10d,'Table10.TE Transition Matrices','A19:F24');


%%       Check the Robustness of Mobility in the Tax Model
%---------------------- Test Different Tax Allowances ---------------------
Mob_dim = 5;          % Dimension of the social mobility matrix (= 5 or 10)
NT_sim  = 301;        % Number of simulations for the test validity
mob_mov = zeros(NT_sim,Mob_dim);  % Mobility probabilities of the quintiles except the top20% & SI
for ss=1:NT_sim
mob_mat = zeros(Mob_dim,Mob_dim,5);   % Social Mobility matrix (mid to end)    
allowance = 0+(ss-1)*0.01;            % Tax allowance relative to income
[sim_W,~,~] = fit_sim(W0,0,0,1,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
sort_W = sort(sim_W);
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=2;
for j=1:Mob_dim % Mobility matrix [row: bottom to top (mid-per), column: bottom to top (end)]
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
for j=1:Mob_dim-1
    mob_mov(ss,j) = 1.0-mob_mat(j,j,k);
end
mob_mov(ss,Mob_dim) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index
clear sim_W;
end

figure(4);
subplot(2,2,1);
   plot(mob_mov(:,1),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.245)*ones(NT_sim,1),'r:','LineWidth',2);             hold on; % baseline mobility
   scatter(101,mob_mov(101,1),90,'o','r','MarkerFaceColor','r');  hold off;
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','1','2','3','4','5','6'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'});
   xlabel(['$Allowance/\overline Y$'],'interpreter','latex','position',[NT_sim-60,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(a) Mobility probability for the initial bottom 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,2);
   plot(mob_mov(:,2),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.240)*ones(NT_sim,1),'r:','LineWidth',2);             hold on;
   scatter(101,mob_mov(101,2),90,'o','r','MarkerFaceColor','r');  hold off;
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','1','2','3','4','5','6'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'});  
   xlabel(['$Allowance/\overline Y$'],'interpreter','latex','position',[NT_sim-60,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(b) Mobility probability for the initial lower 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,3);
   plot(mob_mov(:,3),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.205)*ones(NT_sim,1),'r:','LineWidth',2);             hold on;
   scatter(101,mob_mov(101,3),90,'o','r','MarkerFaceColor','r');  hold off;   
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','1','2','3','4','5','6'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'});   
   xlabel(['$Allowance/\overline Y$'],'interpreter','latex','position',[NT_sim-60,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(c) Mobility probability for the initial middle 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,4);
   plot(mob_mov(:,5),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot(0.955*ones(NT_sim,1),'r:','LineWidth',2);                 hold on;
   scatter(101,mob_mov(101,5),90,'o','r','MarkerFaceColor','r');  hold off;   
   axis([1 NT_sim 0.5 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','1','2','3','4','5','6'});
   xlabel(['$Allowance/\overline Y$'],'interpreter','latex','position',[NT_sim-60,0.45],'Fontname','Times new roman');
   ylabel('Shorrocks Index','Fontname','Times new roman');
   title('(d) Shorrocks Index','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
L=legend({'Tax Model','Baseline Model','$Allowance/\overline Y=1$'},'NumColumns',3,'Position',[0.02,0,1,0.05]);
set(L,'Orientation','horizon','Box','off','FontSize',12,'Fontname','Times new roman','interpreter','latex');
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 12 8]);
print(gcf,'-dpng','-r1080','Figure Tax Mobility Allowance.jpg');

%% -------------------- Test Different Tax Rates ---------------------
allowance = 1;        % Tax allowance relative to income
Mob_dim = 5;          % Dimension of the social mobility matrix (= 5 or 10)
NT_sim = 301;
mob_mov = zeros(NT_sim,Mob_dim);
for ss=1:NT_sim
mob_mat = zeros(Mob_dim,Mob_dim,5);   % Social Mobility matrix (mid to end)    
tax = 0+(ss-1)*0.001;                 % Proportional income tax rate
[sim_W,~,~] = fit_sim(W0,0,0,1,0,0,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance);
sort_W = sort(sim_W);
sort_mid = sortrows(sim_W,mid_per);
sort_mid_nper = sort_mid(:,nper); % End-period wealth sorted by mid-period wealth
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mob_dim-1,1);
for i=1:Mob_dim-1
    bound_btm(i,1) = sort_W(npel/Mob_dim*i,nper);
end
k=2;
for j=1:Mob_dim
    mob_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(1,1))/(npel/Mob_dim);
    for i=2:Mob_dim-1
        mob_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)<=bound_btm(i,1))/(npel/Mob_dim);
    end
    mob_mat(j,Mob_dim,k)=sum(sort_mid_nper((j-1)*npel/Mob_dim+1:j*npel/Mob_dim,1)>bound_btm(Mob_dim-1,1))/(npel/Mob_dim); 
end
for j=1:Mob_dim-1
    mob_mov(ss,j) = 1.0-mob_mat(j,j,k);
end
mob_mov(ss,Mob_dim) = (Mob_dim-trace(mob_mat(:,:,k)))/(Mob_dim-1);  % Static Shorrocks Index
clear sim_W;
end

figure(5);
subplot(2,2,1);
   plot(mob_mov(:,1),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.245)*ones(NT_sim,1),'r:','LineWidth',2);             hold on;
   scatter(201,mob_mov(201,1),90,'o','r','MarkerFaceColor','r');  hold off;
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','10%','20%','30%','40%','50%','60%'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'}); 
   xlabel('Tax rate','position',[NT_sim-30,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(a) Mobility probability for the initial bottom 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,2);
   plot(mob_mov(:,2),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.240)*ones(NT_sim,1),'r:','LineWidth',2);             hold on;
   scatter(201,mob_mov(201,2),90,'o','r','MarkerFaceColor','r');  hold off;
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','10%','20%','30%','40%','50%','60%'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'});  
   xlabel('Tax rate','position',[NT_sim-30,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(b) Mobility probability for the initial lower 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,3);
   plot(mob_mov(:,3),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot((1-0.205)*ones(NT_sim,1),'r:','LineWidth',2);             hold on;
   scatter(201,mob_mov(201,3),90,'o','r','MarkerFaceColor','r');  hold off; 
   axis([1 NT_sim 0 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','10%','20%','30%','40%','50%','60%'});
   set(gca,'ytick',0:0.2:1);
   set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'});    
   xlabel('Tax rate','position',[NT_sim-30,-0.1],'Fontname','Times new roman');
   ylabel('Mobility probability','Fontname','Times new roman');
   title('(c) Mobility probability for the initial middle 20%','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
subplot(2,2,4);
   plot(mob_mov(:,5),'color',[0 0.4470 0.7410],'LineWidth',1.5);  hold on;
   plot(0.955*ones(NT_sim,1),'r:','LineWidth',2);                 hold on;
   scatter(201,mob_mov(201,5),90,'o','r','MarkerFaceColor','r');  hold off;  
   axis([1 NT_sim 0.5 1]);
   set(gca,'xtick',1:100:NT_sim);
   set(gca,'xticklabel',{'0','10%','20%','30%','40%','50%','60%'});
   xlabel('Tax rate','position',[NT_sim-30,0.45],'Fontname','Times new roman');
   ylabel('Shorrocks Index','Fontname','Times new roman');
   title('(d) Shorrocks Index','Fontname','Times new roman');
   set(gca,'FontSize',12,'Fontname','Times new roman');
L=legend({'Tax Model','Baseline Model','Tax Rate = 20%'},'NumColumns',3,'Position',[0.02,0,1,0.05]);
set(L,'Orientation','horizon','Box','off','FontSize',12,'Fontname','Times new roman');
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 12 8]);
print(gcf,'-dpng','-r1080','Figure Tax Mobility Rate.jpg');


%==========================================================================
%%                             Functions
function[WE,ICM]=introduction(W0,Y,nper,C0,C1)
% INPUT--------------------------------------------------------------------
% M1:Labour model, M2:Growth model, M3:Tax model, M4:Capital model
% Mi can take on two values, "1" means being considered while "0" means not
% W0:        Initial wealth at the begining of period 1
% nper:      Number of periods
% c0:        Individual constant consumption level (lower bound)
% c1:        Proportional consumption out of individual income
% OUTPUT-------------------------------------------------------------------
% WE:        Individual Wealth at the end of each period (npel by nper)
% ICM:       Individual total income in each period
npel = size(W0,1);                 % Number of people
WB   = zeros(npel,nper);           % Wealth at the beginning of each period
WE   = zeros(npel,nper);
CSM  = zeros(npel,nper);
ICM  = zeros(npel,nper);
WB(:,1) = W0(:,1);
rng(1);  % Fix the time seed for random generator (this syntax is newest)
for j=1:nper
    y2=zeros(npel,1);      
    % Distribute GDP after taxation: divide into npel shares
    gain0=unidrnd(npel,npel,1); % All the agents have right to have income
    for i=1:npel
        y2(gain0(i,1),1)=y2(gain0(i,1),1)+Y/npel;         
    end    
    ICM(:,j)=y2(:,1); % Total income   
    % Individual consumption "CSM" is   
    CSM(:,j)=C0+C1*(ICM(:,j));
    % Individual wealth at the end of peiod t is
    WE(:,j)=WB(:,j)+ICM(:,j)-CSM(:,j);
    WB(:,j+1)=WE(:,j); 
end
end

%==========================================================================
% Function 2: Simulate data using the Single-factor Model
function[WE,ICM,ICM2]=fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance)
% NEW INPUT----------------------------------------------------------------
% M1:Labour model, M2:Growth model, M3:Tax model, M4:Capital model
% Mi can take on two values, "1" means being considered while "0" means not
% M5:Rent-seeking agents scenario (1% of agents receive privileges)
% W0:        Initial wealth at the begining of period 1
% nper:      Number of periods
% c0:        Individual constant consumption level (lower bound)
% c1:        Proportional consumption out of individual income
% edr:       Mrginal return in terms of human capital
% growth:    A net exogenous growth rate
% OUTPUT ------------------------------------------------------------------
% WE:        Individual Wealth at the end of each period (npel by nper)
% ICM:       Individual total income in each period
npel = size(W0,1);         % npel: Number of people
WB   = zeros(npel,nper+1); % Wealth at the beginning of each period
WE   = zeros(npel,nper);   % Wealth at the end of each period
CSM  = zeros(npel,nper);   % Individual consumption if consumption tax is in
ICM  = zeros(npel,nper);   % Individual total income
ICM2 = zeros(npel,nper);   % Individual income before transfer revenue
WB(:,1)   = W0(:,1);

rent_mask = false(npel,1);
if M5==1
   rent_count = max(1,ceil(0.01*npel));
   rent_mask(1:rent_count,1) = true;
end

rng(1);  % Fix the time seed for random generator (this syntax is newest)
if M1==1
   effort = normrnd(0,sigma,npel,1);
else
   effort = zeros(npel,1);
end
if M2==1
   G = growth;
else
   G = 0.0;
end
for j=1:nper
    needs_subsidy = WB(:,j)<(C0*(1+G)^(j-1));
    num_sub = sum(needs_subsidy); % Number of agents who need subsidy
% Round 1: Subsidy the very poor
    y1 = zeros(npel,1);
    for i=1:npel
        if (M5==1 && rent_mask(i,1)) || WB(i,j) < C0*(1+G)^(j-1)
           y1(i,1) = C0*(1+G)^(j-1);       % The 1st part of income
        end
    end
% Round 2: Distribute the rest output
    if M3==1
       T = tax;
    else
       T = 0;
    end
    if M5==1
       additional_rent = sum(rent_mask & ~needs_subsidy);
    else
       additional_rent = 0;
    end
    rest = Y*(1+G)^(j-1) - (num_sub+additional_rent)*C0*(1+G)^(j-1);
    rest = max(rest,0);
% Creates probabilities determined by Human Capital
    rawprob = zeros(npel,1);
    rawprob(:,1) = 1 - 1./(ones(npel,1)+exp(ones(npel,1)*log(1/(npel-1))+beta1*effort));
    prob1 = rawprob/sum(rawprob); % Logit probability
% Creates probabilities determined by Physical Capital Investment
    inv_share = zeros(npel,1);
    if j==1
       inv_share(:,1) = 1/npel;   % No investment effect in period 1
    else   
       inv_share(:,1) = WE(:,j-1)/sum(WE(:,j-1));
    end   
% Final probability is a weighted sum of Human and Investment
    if M4==0
       prob = prob1;
    elseif M4==1 && M1==0                       
       prob = inv_share;
    elseif M1==1 && M4==1   
       prob = 0.70*prob1 + 0.30*inv_share;
    end    
%------1st, distribute incomeby by npel identical shares    
    gain0 = randsrc(npel,1,[1:npel; prob']); % All the agents have right to have income
    y2 = zeros(npel,1);    
    for i=1:npel
        y2(gain0(i,1),1) = y2(gain0(i,1),1)+rest/npel;         
    end   
%------2nd, charge tax
    taxpay = zeros(npel,1);
    for i=1:npel
        if y2(i,1)==0.0
           taxpay(i,1) = 0.0;
        else
           if y1(i,1)+y2(i,1)>allowance*(Y*(1+G)^(j-1)/npel)
              taxpay(i,1) = (y1(i,1)+y2(i,1)-allowance*(Y*(1+G)^(j-1)/npel))*T;
           else
              taxpay(i,1) = 0.0;
           end
        end
        if M5==1 && rent_mask(i,1)
           taxpay(i,1) = taxpay(i,1)*(1-0.73);
        end
    end
    Taxall = sum(taxpay);
%------3rd, transfer the total tax revenue to agents identically
    ICM(:,j) = y1 + y2 - taxpay + Taxall/npel;     % Total income 
    ICM2(:,j)= y2 - taxpay + Taxall/npel;          % Income before transfer   
%------4th,consider individual consumption "CSM"
    for i=1:npel
        if WB(i,j)+ICM(i,j)-Taxall/npel < C0*(1+G)^(j-1) + C1*(ICM2(i,j)-Taxall/npel)
           CSM(i,j) = WB(i,j) + ICM(i,j)-Taxall/npel;           
        else
           CSM(i,j) = C0*(1+G)^(j-1) + C1*(ICM2(i,j)-Taxall/npel);
        end
    end
% All this revenue will be allocated to agents with same way of the income tax
%------5th, calculate individual wealth after consumption
    WE(:,j)  = WB(:,j) + ICM(:,j) - CSM(:,j);
    WB(:,j+1)= WE(:,j); 
end

end

%==========================================================================
% Function 3: Calculate the Gini coefficients                         
%---Horizontal axi: accumulated sum from 0 to 1 with each interval =1/npel
%   where npel is the number of agents
%---Vertical axi: accumulated wealth shares of population intervals where
%   the total share is already scaled to 1 
function y = gini_coeff(x)
% x has already been the accumlative wealth shares
n = size(x,1);% n is the row number of x
bottom = zeros(n,1);
for i = 2:n 
    bottom(i,1) = (x(i-1,1)+x(i,1))*(1/n)/2;  % calculate the area of a gradiant
end
y = (0.5-sum(bottom(:,1)))/0.5;
end


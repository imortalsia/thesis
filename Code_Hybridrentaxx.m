%%                    Section1: Parameter setting
clear; clc;
% Configuration
npel = 100;       % Number of agents
nper = 100;       % Number of periods
mid_per=nper*0.3;  % The starting period for social mobility 30 years
num_country=8;     % Number of countries [1.Canada;2.China;3.France;4.Germany;5.Italy;6.Japan;7.UK;8.US ]
top = 0.1;         % The percentage for the top share estimation
TM_top = 0.5;      % The percentage (top) for the transition matrix estimation
TM_btm = 0.5;      % The percentage (bottom) for the transition matrix estimation
Mobl_dim = 5;      % Dimension of the social mobility matrix (= 5 or 10)
navg_s = 10;       % Spell Length of Average Gini (short-term): [t,...,t+navg-1]
navg_m = 50;       % Spell Length of Average Gini (mid-term)
navg_l = 100;      % Spell Length of Average Gini (long-term)
save settings_hybrid.mat;  % All settings for later import
[actual_data,~,~] = xlsread('country_data','load_data','A1:L9'); % Read actual data
effort_rates = actual_data(:,1); % Return rate of human capital
growth_rates = actual_data(:,2); % Annual growth rate of GDP
allowances   = actual_data(:,3); % Income tax allowance over GDP per capita
tax_rates    = actual_data(:,4); % Income tax rate for the taxable income
consum0      = actual_data(:,5); % Subsistence consumption level
consum1      = actual_data(:,6); % Consumption propensity
act_giniW    = actual_data(:,7); % Actual Wealth Gini coefficient
act_top10    = actual_data(:,8); % Actual Top 10% wealth share
labor_pwr    = actual_data(:,9); % Labor share of national income
%KY_ratio     = actual_data(:,10);
%labor_pwr    = 0.7*ones(num_country,1);
KY_ratio     = 3.0*ones(num_country,1); % K-Y ratio
tmseed       = 1;                       % Time seed for simulations
optimal_pars = zeros(num_country,1);    % Estimated parameter values
optimal_loss = zeros(num_country,1);    % Loss functions for the estimators

%%                      Section2: Estimation
% Note: You can use "parallel for" if have a high performance PC.
%mypool= parpool(num_country);  
%parfor k=1:num_country 
for k=1:num_country    
    effort    = effort_rates(k,1);
    growth    = growth_rates(k,1);
    tax       = tax_rates(k,1);
    allowance = allowances(k,1);
    C0 = consum0(k,1);   
    C1 = consum1(k,1);
    act_giniw = act_giniW(k,1);    
    act_10    = act_top10(k,1);
    gamma     = labor_pwr(k,1);
    ky_rt     = KY_ratio(k,1);
    tmseed    = 1; % Time seed for simulations
    options = optimset('Display','none','TolX',1e-3,'TolFun',1e-6);
    [x,fval]= fminbnd(@(x) country_optimal(x,C0,C1,effort,growth,tax,allowance,act_giniw,act_10,npel,nper,gamma,ky_rt,tmseed),1,30,options);
    optimal_pars(k,1) = x;
    optimal_loss(k,1) = fval;
end
%delete(mypool);
save optimal_pars.txt optimal_pars -ASCII;
save optimal_loss.txt optimal_loss -ASCII;

%%                    Section 3: Post-estimation
clear; clc;
load settings_hybrid.mat; 
[actual_data,txt,~] = xlsread('country_data','load_data','A1:L9');
sigmas       = load('optimal_pars.txt');
countryname  = txt(2:9,1);
effort_rates = actual_data(:,1);
growth_rates = actual_data(:,2);
allowances   = actual_data(:,3);
tax_rates    = actual_data(:,4);
consum0      = actual_data(:,5);
consum1      = actual_data(:,6);
act_giniW    = actual_data(:,7);
act_top10    = actual_data(:,8);
labor_pwr    = actual_data(:,9);
%KY_ratio     = actual_data(:,10);
%labor_pwr    = 0.7*ones(num_country,1);
KY_ratio     = 3.0*ones(num_country,1);
tmseed = 1; % Time seed for simulations
Y  = npel;       
M1 = 1;     % 1: Consider human capital;    0: Ignore human capital
M2 = 1;     % 1: Consider economic growth;  0: Ignore economic growth
M3 = 1;     % 1: Consider taxation;         0: Ignore taxation
M4 = 1;     % 1: Consider physical capital; 0: Ignore physical capital
M5 = 0;     % 1: Enable rent-seeking agents; 0: Disable rent-seeking agents
if exist('settings_hybrid.mat','file') == 2
    save('settings_hybrid.mat','M5','-append');
end
Lorenz_period= nper;                 % Period at which Lorenz is plotted
wealth_gini  = zeros(num_country,1); % Wealth Gini at "Lorenz_period"
t10s         = zeros(num_country,1); % Top 10% share at "Lorenz_period"
b50s         = zeros(num_country,1); % Bottom 50% share at "Lorenz_period"
up_rate      = zeros(num_country,1); % Upward rate from youth to the old
down_rate    = zeros(num_country,1); % Downward rates from youth to the old
W_sorted_end = zeros(npel,num_country); % Sorted wealth distribution
W_sorted_end1= zeros(npel,num_country); % Scaled sorted wealth distribution 
mobility_mat = zeros(Mobl_dim,Mobl_dim,num_country); % Mobility matrix (mid to end)
Shorrocks_S  = zeros(num_country,1);    % Static Shorrocks Indices
TaxRates     = zeros(nper,num_country); % Tax rate in each period
%mypool= parpool(num_country);
%parfor k=1:num_country
for k=1:num_country
    W0 = ones(npel,1)*KY_ratio(k,1); % Initial wealth endowment
    effort = effort_rates(k,1);
    growth = growth_rates(k,1);
    tax = tax_rates(k,1);
    allowance = allowances(k,1);
    C0 = consum0(k,1);    
    C1 = consum1(k,1);
    sigma = sigmas(k,1);   
    gamma = labor_pwr(k,1);
    gamma_profile = assign_chronic_gamma(gamma, npel, M5);
%-------Call simulation function
    [sim_wealth,~,~,TaxRate] = fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,effort,sigma,growth,tax,allowance,gamma_profile,tmseed);
    sort_wealth   = sort(sim_wealth); 
    TaxRates(:,k) = TaxRate;
    t10s(k,1) = sum(sort_wealth(npel*(1-top)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
    b50s(k,1) = sum(sort_wealth(1:npel*TM_btm,nper),1)/sum(sort_wealth(:,nper),1);    
    W_sorted_end(:,k) = sort_wealth(:,nper);
    W_sorted_end1(:,k)= W_sorted_end(:,k)/mean(W_sorted_end(:,k)); % Scale
%-------Lorenz Curve & Gini
    scale_sortw  = sort_wealth(:,Lorenz_period);                   % Not Scale    
    accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));
    wealth_gini(k,1) = gini_coeff(accum_wealth);       
%-------Calculate upward & downward rates from youth to the old
    Threshold_top = sort_wealth(npel*(1-TM_top),nper);   % Lower bound of the top
    Threshold_btm = sort_wealth(npel*TM_btm+1,nper);     % Upper bound of the bottom
    sort_mid      = sortrows(sim_wealth,mid_per);
    sort_mid_nper = sort_mid(:,nper);  % End-period wealth sorted by the mid-period wealth
    up_rate(k,1)  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
    down_rate(k,1)= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<Threshold_btm)/(npel*TM_top);    
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
    bound_btm = zeros(Mobl_dim-1,1);
    for i=1:Mobl_dim-1
        bound_btm(i,1) = sort_wealth(npel/Mobl_dim*i,nper);
    end
    % Mobility matrix [row: bottom to top at the beginning, column: bottom to top at then end]
    for j=1:Mobl_dim
        mobility_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(1,1))/(npel/Mobl_dim);
        for i=2:Mobl_dim-1
            mobility_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(i-1,1) ...
            & sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(i,1))/(npel/Mobl_dim);
        end
        mobility_mat(j,Mobl_dim,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(Mobl_dim-1,1))/(npel/Mobl_dim); 
    end
    Shorrocks_S(k,1) = (Mobl_dim-trace(mobility_mat(:,:,k)))/(Mobl_dim-1);
end
%delete(mypool);

%%               Make Table 6 & 11
countries  = {'Canada';'China';'France';'Germany';'Italy';'Japan';'UK';'US'};
table1     = table(countries, wealth_gini, t10s, up_rate, Shorrocks_S, ...
    'VariableNames', {'Country','Gini','Top10','Up','SI'});

filename_tables = 'Generated Tables.xlsx';
sheet_table7    = 'Table7.Inequality and Mobility';
sheet_table11   = 'Table11.Mobility Matrices';

if isfile(filename_tables)
    existingSheets = {};
    if exist('sheetnames','file') == 2 || exist('sheetnames','builtin') == 5
        try
            existingSheets = sheetnames(filename_tables);
        catch
            existingSheets = {};
        end
    end
    if ~any(strcmp(sheet_table7, existingSheets))
        delete(filename_tables);
    end
end

writeArgs = {'Sheet', sheet_table7};
if isfile(filename_tables)
    writeArgs = [writeArgs, {'WriteMode','overwritesheet'}];
end

try
    writetable(table1, filename_tables, writeArgs{:});
catch ME
    if contains(ME.message, 'Unrecognized parameter name') || contains(ME.identifier, 'UnrecognizedParameter')
        writetable(table1, filename_tables, 'Sheet', sheet_table7);
    else
        rethrow(ME);
    end
end

table2Bcolm = {'Bottom 20%';'Lower 20%';'Middle 20%';'Upper 20%';'Top 20%'};
table2Bvals = num2cell(Mobl_dim,Mobl_dim);
Ncolm2B     = Mobl_dim+1;
table2Brows = num2cell(1,Ncolm2B);
table_app2B = num2cell(num_country*Ncolm2B,Ncolm2B);
for i=1:num_country
    table2Bvals = num2cell(mobility_mat(:,:,i));
    table2Brows = [countries(i,1),table2Bcolm'];
    table_app2B((i-1)*Ncolm2B+1,1:Ncolm2B) = table2Brows;
    table_app2B((i-1)*Ncolm2B+2:(i-1)*Ncolm2B+Ncolm2B,1:Ncolm2B) = [table2Bcolm,table2Bvals];
end

writeCellArgs = {'Sheet', sheet_table11, 'WriteMode','overwritesheet'};
try
    writecell(table_app2B, filename_tables, writeCellArgs{:});
catch ME
    if contains(ME.message, 'Unrecognized parameter name') || contains(ME.identifier, 'UnrecognizedParameter')
        writecell(table_app2B, filename_tables, 'Sheet', sheet_table11);
    else
        rethrow(ME);
    end
end
%%        Figure 4.Simulated and Observed Inequality
figHandle = create_new_figure('Figure4.Simulated and Observed Wealth Inequality');
   scat1 = scatter(act_giniW,act_top10,50,'o','k','MarkerFaceColor','k');
   for i=1:num_country
	   text(act_giniW(i)-0.011,act_top10(i)+0.004,countryname(i),'FontName','Times','Fontsize',12,'Color','k');
       text(wealth_gini(i)-0.002,t10s(i)-0.009,countryname(i),'FontName','Times','Fontsize',12,'Color','b');
   end
   hold on;
   scat2 = scatter(wealth_gini,t10s,60,'d','b','MarkerFaceColor','b');
   clear i;  % Now plot the linking line between nodes
   for i=1:num_country
       hold on;
       nodex = [act_giniW(i); wealth_gini(i)];
       nodey = [act_top10(i); t10s(i)];
       plot(nodex,nodey,':r','LineWidth',2);
   end  
   hold off;
   axis([0.6 0.9 0.4 0.8]);
   xlabel('Wealth Gini','FontSize',12,'FontName','Times New Roman');
   ylabel('Top10% Wealth Share','FontSize',12,'FontName','Times New Roman');  
   legend([scat1,scat2],'Acutal','Simulated','Location','NorthWest','FontName','Times New Roman','Fontsize',12);
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 9 8]);
print(figHandle,'-dpng','-r1080','Figure4.Simulated and Observed Wealth Inequality.jpg');

%%        Figure 7.Wealth Distributions of the Rich & the Rest (US)
country_test = 8;                              % Country to be tested (US)
ntop_test  = npel*0.2;                         % Whom to be tested
top_test_x = (ntop_test:-1:1)';                % Name the order of the rich
top_test_y = W_sorted_end1(npel-ntop_test+1:npel,country_test); 
btm_test_x = (npel-ntop_test:-1:+1)';          % Name the order of the poor
btm_test_y = W_sorted_end1(1:npel-ntop_test,country_test); 
[xData1,yData1] = prepareCurveData(top_test_x,top_test_y);
[xData2,yData2] = prepareCurveData(btm_test_x,btm_test_y);
ft_power = fittype('power1');                              % Set up fittype
ft_exp   = fittype( 'exp1' );
opts = fitoptions('Method','NonlinearLeastSquares');       % Set up options
opts.Display = 'Off';
figHandle = create_new_figure('Figure7.Wealth Distributions of the Rich and the Rest (US)'); % Plot the fitted curves of simulated distribution
subplot(2,2,1)
    opts.StartPoint = [10000 -0.7];
    [fitresult1,gof1] = fit(xData1,yData1,ft_power,opts);  % Fit to data
    hfit1 = plot(fitresult1,xData1,yData1);
    set(hfit1,'LineWidth',2,'MarkerSize',12);
    axis_value=axis;
    axis([axis_value(1) axis_value(2) 0 60]);
    xlabel('Rank of individual wealth','FontName','Times New Roman');
    ylabel('Standardized individual wealth','FontName','Times New Roman');
    legend off;
    interval=get(gca);
    x_itv=interval.XLim;              % obtain upper & lower limit of x axis
    y_itv=interval.YLim;              % obtain upper & lower limit of y axis
    pot=[0.35 0.85];                  % set up relative position of the text
    x_itv0=x_itv(1)+pot(1)*(x_itv(2)-x_itv(1));   % horizontal value of text
    y_itv0=y_itv(1)+pot(2)*(y_itv(2)-y_itv(1));   % vertical value of text
% Automatically read estimaters and write functions 
    coeffit = string(coeffvalues(fitresult1));    
    coef1=sprintf('%5.2f',coeffit(1));   coef2=sprintf('%5.3f',coeffit(2));   
    tstr1 = {'\itY  = '};  tstr2 = strrep({'*\itX^{power}'},'power',coef2);
    tstr3 = strrep({'R^2= rsq'},'rsq',sprintf('%5.3f',gof1.rsquare));
    tstr = [strcat(tstr1,coef1,tstr2);tstr3];
    text(x_itv0,y_itv0,tstr,'edgecolor','w','FontSize',12);
    title('(a1) Power Law fitting for the rich','FontName','Times New Roman');
    set(gca,'Fontsize',12);
subplot(2,2,2)
    opts.StartPoint = [3600 -0.015];
    [fitresult2,gof2] = fit(xData1,yData1,ft_exp,opts);       % Fit to data
    hfit2 = plot(fitresult2,xData1,yData1);
    set(hfit2,'LineWidth',2,'MarkerSize',12);
    axis_value=axis;
    axis([axis_value(1) axis_value(2) 0 60]);
    xlabel('Rank of individual wealth','FontName','Times New Roman');
    ylabel('Standardized individual wealth','FontName','Times New Roman');
    legend off;
    interval=get(gca);
    x_itv=interval.XLim;              % obtain upper & lower limit of x axis
    y_itv=interval.YLim;              % obtain upper & lower limit of y axis
    pot=[0.35 0.85];                  % set up relative position of the text
    x_itv0=x_itv(1)+pot(1)*(x_itv(2)-x_itv(1));   % horizontal value of text
    y_itv0=y_itv(1)+pot(2)*(y_itv(2)-y_itv(1));   % vertical value of text
    coeffit = string(coeffvalues(fitresult2));
    coef1=sprintf('%5.2f',coeffit(1));   coef2=sprintf('%5.3f',coeffit(2));   
    tstr1={'\itY  = '};  tstr2=strrep({'*exp({power}\itX)'},'power',coef2);
    tstr3 = strrep({'R^2= rsq'},'rsq',sprintf('%5.3f',gof2.rsquare));
    tstr = [strcat(tstr1,coef1,tstr2);tstr3];
    text(x_itv0,y_itv0,tstr,'edgecolor','w','FontSize',12);
    title('(a2) Exponential fitting for the rich','FontName','Times New Roman');
    set(gca,'Fontsize',12); 
subplot(2,2,3)
    opts.StartPoint = [1045 -0.008];
    [fitresult3,gof3] = fit(xData2,yData2,ft_exp,opts);      % Fit to data
    hfit3 = plot(fitresult3,xData2,yData2);
    set(hfit3,'LineWidth',2,'MarkerSize',12);
    axis_value=axis;
    axis([axis_value(1) axis_value(2) 0 1]);
    xlabel('Rank of individual wealth','FontName','Times New Roman');
    ylabel('Standardized individual wealth','FontName','Times New Roman');
    legend off;
    interval=get(gca);
    x_itv=interval.XLim;              % obtain upper & lower limit of x axis
    y_itv=interval.YLim;              % obtain upper & lower limit of y axis
    pot=[0.35 0.85];                  % set up relative position of the text
    x_itv0=x_itv(1)+pot(1)*(x_itv(2)-x_itv(1));   % horizontal value of text
    y_itv0=y_itv(1)+pot(2)*(y_itv(2)-y_itv(1));   % vertical value of text
    coeffit = string(coeffvalues(fitresult3));
    coef1=sprintf('%5.3f',coeffit(1));   coef2=sprintf('%5.3f',coeffit(2));   
    tstr1={'\itY  = '};  tstr2=strrep({'*exp({power}\itX)'},'power',coef2);
    tstr3 = strrep({'R^2= rsq'},'rsq',sprintf('%5.3f',gof3.rsquare));
    tstr = [strcat(tstr1,coef1,tstr2);tstr3];
    text(x_itv0,y_itv0,tstr,'edgecolor','w','FontSize',12);
    title('(b1) Exponential fitting for the poor','FontName','Times New Roman');
    set(gca,'Fontsize',12);
subplot(2,2,4)
    opts.StartPoint = [1044089875 -3];
    [fitresult4,gof4] = fit(xData2,yData2,ft_power,opts);     % Fit to data
    hfit4 = plot(fitresult4,xData2,yData2);
    set(hfit4,'LineWidth',2,'MarkerSize',12);
    axis_value=axis;
    axis([axis_value(1) axis_value(2) 0 1]);
    xlabel('Rank of individual wealth','FontName','Times New Roman');
    ylabel('Standardized individual wealth','FontName','Times New Roman');
    legend off;
    interval=get(gca);
    x_itv=interval.XLim;              % obtain upper & lower limit of x axis
    y_itv=interval.YLim;              % obtain upper & lower limit of y axis
    pot=[0.35 0.85];                  % set up relative position of the text
    x_itv0=x_itv(1)+pot(1)*(x_itv(2)-x_itv(1));   % horizontal value of text
    y_itv0=y_itv(1)+pot(2)*(y_itv(2)-y_itv(1));   % vertical value of text
    coeffit = string(coeffvalues(fitresult4));
    coef1=sprintf('%5.3f',coeffit(1));   coef2=sprintf('%5.3f',coeffit(2));   
    tstr1 = {'\itY  = '};  tstr2 = strrep({'*\itX^{power}'},'power',coef2);
    tstr3 = strrep({'R^2= rsq'},'rsq',sprintf('%5.3f',gof4.rsquare));
    tstr = [strcat(tstr1,coef1,tstr2);tstr3];
    text(x_itv0,y_itv0,tstr,'edgecolor','w','FontSize',12);    
    title('(b2) Power Law fitting for the poor','FontName','Times New Roman');
    set(gca,'Fontsize',12);   
L = legend(hfit4,'Simulated wealth distribution ',' Fitted curve');
    L.Position = [0.455 0.005 0.1 0.02]; set(L,'Orientation','horizontal','Box','off','FontSize',12,'FontName','Times New Roman');
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 10 8]);
print(figHandle,'-dpng','-r1080','Figure7.Wealth Distributions of the Rich and the Rest (US).jpg');

%%        Figure 8 & Table 6. Estimated Thickness Indices
index_start = npel*0.01+1;             % Top range (ie. top1%) started from
index_end   = npel*0.2+1;              % Top range (ie. top20%) ended up
index_every = 1;                       % Agents added in each rolling 
allreg=size((index_start:index_every:index_end),2); % Regressions number
regB  = zeros(2,allreg,num_country);                % Regression coefs
RegSt = zeros(4,allreg,num_country);                % Regression statistics
clear k;
W_sorted_reverse = sort(W_sorted_end1,'descend');   % Reverse sorting
for nreg = 1:allreg
    clear regx regy RegSt0;   
    regy = log((1:index_start-1+(nreg-1)*index_every)');        
    regx(:,1) = ones(index_start-1+(nreg-1)*index_every,1);
    for k = 1:num_country
        regx(:,2) = log(W_sorted_reverse(1:index_start-1+(nreg-1)*index_every,k));
        [regB(:,nreg,k),~,~,~,RegSt0] = regress(regy,regx);
        format short g;RegSt(:,nreg,k) = (roundn(RegSt0,-4))';
    end        
end
figHandle = create_new_figure('Figure8.Estimates of The Tail Thickness Indices');
    plot((1:allreg),-1*regB(2,:,1),':k','LineWidth',2.5); hold on;
    plot((1:allreg),-1*regB(2,:,2),':b','LineWidth',2.5); hold on;
    plot((1:allreg),-1*regB(2,:,3),':c','LineWidth',2.5); hold on;
    plot((1:allreg),-1*regB(2,:,4),':r','LineWidth',2.5); hold on;    
    plot((1:allreg),-1*regB(2,:,5),'-k','LineWidth',1.5); hold on;    
    plot((1:allreg),-1*regB(2,:,6),'-b','LineWidth',1.5); hold on;
    plot((1:allreg),-1*regB(2,:,7),'-c','LineWidth',1.5); hold on;    
    plot((1:allreg),-1*regB(2,:,8),'-r','LineWidth',1.5); hold off;
    grid on;
    axis_value=axis;
    axis([1 allreg axis_value(3) axis_value(4)]);
    set(gca,'xtick',[1,npel*0.04+1,npel*0.09+1,npel*0.14+1,npel*0.19+1]);
    set(gca,'xticklabel',{'Top1%','Top5%','Top10%','Top15%','Top20%'},'FontName','Times New Roman');
    xlabel('Sample used for estimation','FontName','Times New Roman');
    ylabel('Estimated Thickness Index {\bf\phi}','FontName','Times New Roman');
    set(gca,'Fontsize',12); 
L = legend('Canada ','China ','France ','Germany ','Italy ','Japan ','UK ','US');
    L.Position = [0.460 0.005 0.1 0.025]; set(L,'Orientation','horizontal','Box','off','FontSize',12,'FontName','Times New Roman');
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 10 7]);
print(figHandle,'-dpng','-r1080','Figure8.Estimates of The Tail Thickness Indices.jpg');

% Table6.Estimated Thickness Indices
Index_all = zeros(num_country,4);
clear i;
for i=1:num_country
    Index_all(i,1) =-1*regB(2,1,i);                   % Top 1%
    Index_all(i,2) =-1*regB(2,round(npel*0.02+1),i);  % Top 3%
    Index_all(i,3) =-1*regB(2,round(npel*0.09+1),i);  % Top 10%
    Index_all(i,4) =-1*regB(2,round(npel*0.19+1),i);  % Top 20%
end
table3_header = {'Country','Top1%','Top3%','Top10%','Top20%'};
table3 = table(string(countries(:)),Index_all(:,1),Index_all(:,2),Index_all(:,3),Index_all(:,4), ...
    'VariableNames',matlab.lang.makeValidName(table3_header));
writecell(table3_header,'Generated Tables.xlsx','Sheet','Table6.Estimated Thickness','Range','A1:E1');
writetable(table3,'Generated Tables.xlsx','Sheet','Table6.Estimated Thickness','Range','A2','WriteVariableNames',false);

%%        Figure 9. Intra-generational Great Gatsby Curve
Corr_Gatsby = corr(wealth_gini,Shorrocks_S);
OLS_x = wealth_gini;  
OLS_y = Shorrocks_S;
mdl1  = fitlm(table(OLS_x,OLS_y),'Intercept',true);   % Regression: Yt=β0+β1*Xt+e
OLS_CoefTab = mdl1.Coefficients;   OLS_Coef = OLS_CoefTab.Variables;
B_value1    = OLS_Coef(:,1);  
SE_value1   = OLS_Coef(:,2);  
P_value1    = OLS_Coef(:,4);
xData5_ext  = (0.6:0.1:0.9)';
yData5_ext  = B_value1(1,1)+B_value1(2,1)*xData5_ext;
figHandle = create_new_figure('Figure9.Intra-generational Great Gatsby Curve');
   scat3 = scatter(wealth_gini,Shorrocks_S,80,'filled','b','MarkerFaceColor','b'); hold on;
   hfit5 = plot(xData5_ext,yData5_ext,'r','LineWidth',2); hold off;
   for i=1:num_country
       text(wealth_gini(i)-0.008,Shorrocks_S(i)-0.016,countryname(i),'FontName','Times','Fontsize',12,'Color','k');
   end
   interval=get(gca);
   x_itv=interval.XLim;              % obtain upper & lower limit of x axis
   y_itv=interval.YLim;              % obtain upper & lower limit of y axis
   pot=[0.9 0.95];                   % set up relative position of the text
   x_itv0=x_itv(1)+pot(1)*(x_itv(2)-x_itv(1));   % horizontal value of text
   y_itv0=y_itv(1)+pot(2)*(y_itv(2)-y_itv(1));   % vertical value of text
   coefstr = string(OLS_Coef);
   coef1=sprintf('%5.3f',coefstr(1,1));   coef2=sprintf('%5.3f',coefstr(2,1));   
   tstr1={'\itY = beta1  + beta2*\itX'};  
   tstr2=strrep(tstr1,'beta1',coef1); tstr3=strrep(tstr2,'beta2',coef2);
   tstr4 = strrep({'       (SE1)  (SE2)'},'SE1',sprintf('%5.3f',coefstr(1,2)));
   tstr5 = strrep(tstr4,'SE2',sprintf('%5.3f',coefstr(2,2)));
   tstr = [tstr3;tstr5];
   text(x_itv0,y_itv0,tstr,'edgecolor','w','FontSize',12);
   axis([0.6 1.0 0.2 0.7]);
   xticks(0.6:0.1:1);     xticklabels(0.6:0.1:1);
   yticks(0.2:0.1:0.7);   yticklabels(0.2:0.1:0.7);
   xlabel('Wealth Gini Coefficients','FontSize',12,'FontName','Times New Roman');
   ylabel('Shorrocks Index','FontSize',12,'FontName','Times New Roman'); 
   legend([scat3,hfit5],'Simulated Data','Fitted Line','Orientation','vertical','Location','Northwest','Box','off','FontSize',12);
   set(gca,'Fontsize',12,'FontName','Times New Roman');
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 9 6.5]);
print(figHandle,'-dpng','-r1080','Figure9.Intra-generational Great Gatsby Curve.jpg');


%%        Robustness
% Reset the "npel" (Number of agents) or "nper" (Number of periods) in
% section 1 and rerun the steps above.


%%        Monte Carlo 1: Distribution of simulated Gini & Top
% Simulate distributions of Gini & Top shares based on the estimated sigma
% Steps: Directly run this section
clear; clc;
load settings_hybrid.mat;
[actual_data,txt,~] = xlsread('country_data','load_data','A1:L9');
countryname  = txt(2:9,1);
effort_rates = actual_data(:,1);
growth_rates = actual_data(:,2);
allowances   = actual_data(:,3);
tax_rates    = actual_data(:,4);
consum0      = actual_data(:,5);
consum1      = actual_data(:,6);
act_giniW    = actual_data(:,7);
act_top10    = actual_data(:,8);
labor_pwr    = actual_data(:,9);
%KY_ratio     = actual_data(:,10);
%labor_pwr    = 0.7*ones(num_country,1);
KY_ratio     = 3.0*ones(num_country,1);
sigmas       = load('optimal_pars.txt');
Y  = npel;       
M1 = 1;
M2 = 1;
M3 = 1;
M4 = 1;
M5 = 1;
Lorenz_period = nper;
MCnum = 100;                                    % Number of MC simulations
MC_gini_all = zeros(MCnum,num_country);          % Save MC parameter values
MC_top_all  = zeros(MCnum,num_country); 
% Run MCnum sample groups
for smn=1:MCnum
    tmseed       = smn;
    wealth_gini  = zeros(num_country,1);
    t10s         = zeros(num_country,1);
    b50s         = zeros(num_country,1);
    W_sorted_end = zeros(npel,num_country);
    W_sorted_end1= zeros(npel,num_country);
%    mypool= parpool(num_country);
%    parfor k=1:num_country
    for k=1:num_country
        W0     = ones(npel,1)*KY_ratio(k,1);  
        effort = effort_rates(k,1);
        growth = growth_rates(k,1);
        tax    = tax_rates(k,1);
        allowance = allowances(k,1);
        C0    = consum0(k,1);    
        C1    = consum1(k,1);
        sigma = sigmas(k,1);
        gamma = labor_pwr(k,1);
        gamma_profile = assign_chronic_gamma(gamma, npel, M5);
%-------Call simulation function
        [sim_wealth,~,~,~] = fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,effort,sigma,growth,tax,allowance,gamma_profile,tmseed);
        sort_wealth = sort(sim_wealth);   
        t10s(k,1)   = sum(sort_wealth(npel*(1-top)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
        b50s(k,1)   = sum(sort_wealth(1:npel*TM_btm,nper),1)/sum(sort_wealth(:,nper),1);    
        W_sorted_end(:,k) = sort_wealth(:,nper);
        W_sorted_end1(:,k)= W_sorted_end(:,k)/mean(W_sorted_end(:,k));
%-------Lorenz Curve & Gini
        scale_sortw  = sort_wealth(:,Lorenz_period);  
        accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));
        wealth_gini(k,1) = gini_coeff(accum_wealth);  
        
        MC_gini_all(smn,k) = wealth_gini(k,1);
        MC_top_all(smn,k)  = t10s(k,1);
    end
%delete(mypool);
end
save MC_gini_all.txt MC_gini_all -ASCII;
save MC_top_all.txt MC_top_all -ASCII;

%------ Load simulated samples of Gini and Top10, and the actual data
MC_gini_all = load('MC_gini_all.txt');
MC_top_all  = load('MC_top_all.txt');
Nsample     = size(MC_gini_all,1);
%------ Test if simulated Gini & Top10 follow beta distributions
PD_gini  = zeros(num_country,2);
PD_top   = zeros(num_country,2);
CDF_gini = zeros(num_country,1);   % CDF value at the spot of actual gini
CDF_top  = zeros(num_country,1);   % CDF value at the spot of actual top10
Prb_gini = zeros(num_country,1);   % Probability of actual gini falling-in
Prb_top  = zeros(num_country,1);   % Probability of actual top10 falling-in
for i=1:num_country
    [PD_gini(i,:),~] = betafit(MC_gini_all(:,i),0.01);
    [PD_top(i,:),~]  = betafit(MC_top_all(:,i),0.01);
    CDF_gini(i,1)    = betacdf(act_giniW(i,1),PD_gini(i,1),PD_gini(i,2));
    if CDF_gini(i,1) > 0.5
       Prb_gini(i,1) = (1.0-CDF_gini(i,1))*2;
    else
       Prb_gini(i,1) = CDF_gini(i,1)*2;
    end
    CDF_top(i,1)    = betacdf(act_top10(i,1),PD_top(i,1),PD_top(i,2));
    if CDF_top(i,1) > 0.5
       Prb_top(i,1) = (1.0-CDF_top(i,1))*2;     
    else
       Prb_top(i,1) = CDF_top(i,1)*2;       
    end    
end
tableMC1_header = {'Country','P-value for Gini','P-value for Top10%'};
tableMC1 = table(countryname, Prb_gini, Prb_top);
writecell(tableMC1_header,'Generated Tables.xlsx','Sheet','Table5.P-Values of Gini & Top','Range','A1:C1');
writetable(tableMC1,'Generated Tables.xlsx','Sheet','Table5.P-Values of Gini & Top','Range','A2','WriteVariableNames',false);

%------ Plot distributions of Gini
country_plot1 = 8; % Country to plot (US)
  ginix1_min = min(MC_gini_all(:,country_plot1));
  ginix1_max = max(MC_gini_all(:,country_plot1));
  ginix1_dif = (ginix1_max-ginix1_min)/(Nsample-1);
  fit1_x = [ginix1_min:ginix1_dif:ginix1_max]';
  fit1_y = betapdf(fit1_x,PD_gini(country_plot1,1),PD_gini(country_plot1,2));
  CRT_L1 = betainv(0.025,PD_gini(country_plot1,1),PD_gini(country_plot1,2));
  CRT_R1 = betainv(0.975,PD_gini(country_plot1,1),PD_gini(country_plot1,2));
country_plot2 = 2; % Country to plot (China)
  ginix2_min = min(MC_gini_all(:,country_plot2));
  ginix2_max = max(MC_gini_all(:,country_plot2));
  ginix2_dif = (ginix2_max-ginix2_min)/(Nsample-1);
  fit2_x = [ginix2_min:ginix2_dif:ginix2_max]';
  fit2_y = betapdf(fit2_x,PD_gini(country_plot2,1),PD_gini(country_plot2,2));
  CRT_L2 = betainv(0.025,PD_gini(country_plot2,1),PD_gini(country_plot2,2));
  CRT_R2 = betainv(0.975,PD_gini(country_plot2,1),PD_gini(country_plot2,2));
ngrid      = 20;
[fp1,xp1]  = ecdf(MC_gini_all(:,country_plot1));
  his1_x   = min(xp1):(max(xp1)-min(xp1))/(ngrid-1):max(xp1);
  his1_y   = ecdfhist(fp1,xp1,ngrid);
  barcolr1 = [.3 .75 .93];
[fp2,xp2]  = ecdf(MC_gini_all(:,country_plot2));
  his2_x   = min(xp2):(max(xp2)-min(xp2))/(ngrid-1):max(xp2);
  his2_y   = ecdfhist(fp2,xp2,ngrid);
  barcolr2 = [1 .4 .4];  
figHandle = create_new_figure('Figure6.MC Distribution of Gini');
h1a=bar(his1_x,his1_y,'FaceColor',barcolr1,'EdgeColor',barcolr1,'LineWidth',1,'FaceAlpha',0.5);
  hold on;
  h1b=plot(fit1_x,fit1_y,'color','blue','LineWidth',3);
  hold on;
  dash_yL = betapdf(CRT_L1,PD_gini(country_plot1,1),PD_gini(country_plot1,2));
  plot(CRT_L1*ones(101,1),[0:dash_yL/100:dash_yL]','b--','LineWidth',3);  
  hold on;
  dash_yR = betapdf(CRT_R1,PD_gini(country_plot1,1),PD_gini(country_plot1,2));
  plot(CRT_R1*ones(101,1),[0:dash_yR/100:dash_yR]','b--','LineWidth',3);  
  hold on;  
  h1c=scatter(act_giniW(country_plot1,1),betapdf(act_giniW(country_plot1,1),PD_gini(country_plot1,1),PD_gini(country_plot1,2)),400,'o','b','MarkerFaceColor','b');
  hold on;
h2a=bar(his2_x,his2_y,'FaceColor',barcolr2,'EdgeColor',barcolr2,'LineWidth',1,'FaceAlpha',0.5);
  hold on;
  h2b=plot(fit2_x,fit2_y,'color','red','LineWidth',3);
  hold on;
  dash_yL = betapdf(CRT_L2,PD_gini(country_plot2,1),PD_gini(country_plot2,2));
  plot(CRT_L2*ones(101,1),[0:dash_yL/100:dash_yL]','r--','LineWidth',3);
  hold on;
  dash_yR = betapdf(CRT_R2,PD_gini(country_plot2,1),PD_gini(country_plot2,2));
  plot(CRT_R2*ones(101,1),[0:dash_yR/100:dash_yR]','r--','LineWidth',3);
  hold on;
  h2c=scatter(act_giniW(country_plot2,1),betapdf(act_giniW(country_plot2,1),PD_gini(country_plot2,1),PD_gini(country_plot2,2)),400,'o','r','MarkerFaceColor','r');
  hold off; 
axis_value=axis;
  axis([axis_value(1) axis_value(2) 0 max(max(his1_y),max(his2_y))*1.1]);
  xlabel('Simulated Gini');
  ylabel('Frequency');
text(CRT_L1*0.990,max(fit1_y)*0.25,'CDF=2.5%','FontName','Times','Fontsize',13,'Color','b','FontWeight','bold');
  text(CRT_R1*0.998,max(fit1_y)*0.25,'CDF=97.5%','FontName','Times','Fontsize',13,'Color','b','FontWeight','bold');
  text(CRT_L2*0.989,max(fit1_y)*0.25,'CDF=2.5%','FontName','Times','Fontsize',13,'Color','r','FontWeight','bold');
  text(CRT_R2*0.997,max(fit1_y)*0.25,'CDF=97.5%','FontName','Times','Fontsize',13,'Color','r','FontWeight','bold');
set(gca,'FontSize',12,'FontName','Times New Roman');
  L=legend([h1a,h1b,h1c,h2a,h2b,h2c],'Histogram of simulated Gini (US)','Fitted Beta pdf (US)','Actual Gini (US)','Histogram of simulated Gini (China)','Fitted Beta pdf (China)','Actual Gini (China)','location','southoutside','NumColumns',3);
  set(L,'Orientation','horizontal','Box','off','FontSize',12,'FontName','Times New Roman');
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 12 7]);
print(figHandle,'-dpng','-r1080','Figure6.MC Distribution of Gini.jpg');

%%        Monte Carlo 2: Distribution of simulated sigmas 
% Treat simulated Gini & Top shares as the true values and reestimate sigma
% Steps: Directly run this section
clear;  clc;
load settings_hybrid.mat;
[actual_data,txt,~] = xlsread('country_data','load_data','A1:L9');
sigmas       = load('optimal_pars.txt');
effort_rates = actual_data(:,1);
growth_rates = actual_data(:,2);
allowances   = actual_data(:,3);
tax_rates    = actual_data(:,4);
consum0      = actual_data(:,5);
consum1      = actual_data(:,6);
act_giniW    = actual_data(:,7);
act_top10    = actual_data(:,8);
labor_pwr    = actual_data(:,9);
%KY_ratio     = actual_data(:,10);
%labor_pwr    = 0.7*ones(num_country,1);
KY_ratio     = 3.0*ones(num_country,1);
ad_sheet_name = 'Table7.Inequality and Mobility';
ad_range      = 'A1:C9';
script_dir    = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
ad_workbook = fullfile(script_dir,'Generated Tables.xlsx');
if exist(ad_workbook,'file') ~= 2
    error('Code_Hybrid_Model:MissingWorkbook', ...
        'Required workbook %s was not found.', ad_workbook);
end

[~, workbook_sheets] = xlsfinfo(ad_workbook);
sheet_candidates = {ad_sheet_name};
if ~isempty(workbook_sheets)
    trimmed_names = strtrim(workbook_sheets);
    sheet_match = strcmpi(trimmed_names, strtrim(ad_sheet_name));

    if ~any(sheet_match)
        % Apply a more permissive matching rule that ignores punctuation and
        % repeated whitespace. Some MATLAB/Octave versions appear to pad sheet
        % names with non-breaking spaces or replace punctuation, which prevents
        % direct string comparisons from succeeding even though the intended
        % sheet exists.
        normalized_requested = normalize_sheet_name(ad_sheet_name);
        normalized_available = cellfun(@normalize_sheet_name, workbook_sheets, ...
            'UniformOutput', false);
        sheet_match = strcmp(normalized_available, normalized_requested);

        if ~any(sheet_match)
            % As a final heuristic, allow partial matches (e.g., when the sheet
            % name is truncated in older Excel readers). Prioritize the first
            % match reported by xlsfinfo to preserve deterministic behaviour.
            partial_match = find(contains(normalized_available, normalized_requested) | ...
                contains(normalized_requested, normalized_available)); %#ok<FNDSB>
            if ~isempty(partial_match)
                sheet_match(partial_match(1)) = true;
            end
        end
    end

    if any(sheet_match)
        sheet_candidates = workbook_sheets(find(sheet_match,1));
    elseif numel(workbook_sheets) >= 2
        sheet_candidates{end+1} = workbook_sheets{2}; %#ok<AGROW>
    end
else
    sheet_candidates{end+1} = 2; %#ok<AGROW>
end

ad_sheet_to_read = [];
read_error = [];
for candidate_idx = 1:numel(sheet_candidates)
    candidate = sheet_candidates{candidate_idx};
    if ischar(candidate)
        candidate_to_try = strtrim(candidate);
    else
        candidate_to_try = candidate;
    end
    try
        AD_true_data = read_sheet_with_fallback(ad_workbook, candidate_to_try, ad_range);
        ad_sheet_to_read = candidate_to_try;
        read_error = [];
        break
    catch read_error
        AD_true_data = [];
    end
end

if isempty(ad_sheet_to_read)
    if isempty(read_error)
        error('Code_Hybrid_Model:MissingWorksheet', ...
            ['Worksheet "%s" was not found in %s, and no suitable fallback ', ...
            'sheet is available.'], ad_sheet_name, ad_workbook);
    else
        rethrow(read_error);
    end
elseif ischar(ad_sheet_to_read) && ~strcmpi(ad_sheet_to_read, ad_sheet_name)
    warning('Code_Hybrid_Model:UsingFallbackSheet', ...
        'Worksheet "%s" was not found. Using "%s" instead.', ...
        ad_sheet_name, ad_sheet_to_read);
elseif isnumeric(ad_sheet_to_read) && ad_sheet_to_read ~= 1
    warning('Code_Hybrid_Model:UsingFallbackSheet', ...
        'Worksheet "%s" was not found. Using sheet index %d instead.', ...
        ad_sheet_name, ad_sheet_to_read);
end

AD_true_gini = AD_true_data(:,1);
AD_true_top  = AD_true_data(:,2);
MCnum = 100;                                % Number of MC simulations
MC_pars_all  = zeros(MCnum,num_country);     % Save the MC parameter values
MC_loss_all  = zeros(MCnum,num_country); 
% Estimation
for smn=1:MCnum
    tmseed  = smn;
    MC_pars = zeros(1,num_country);
    MC_loss = zeros(1,num_country);
%mypool= parpool(num_country);
%parfor k=1:num_country
    for k=1:num_country    
        effort = effort_rates(k,1);  
        growth = growth_rates(k,1);
        tax    = tax_rates(k,1);
        allowance = allowances(k,1);
        C0 = consum0(k,1);  
        C1 = consum1(k,1);
        act_giniw = AD_true_gini(k,1);    
        act_10    = AD_true_top(k,1);
        gamma = labor_pwr(k,1);
        ky_rt = KY_ratio(k,1);    
        options  = optimset('Display','none','TolX',1e-3,'TolFun',1e-6);
        [x,fval] = fminbnd(@(x) country_optimal(x,C0,C1,effort,growth,tax,allowance,act_giniw,act_10,npel,nper,gamma,ky_rt,tmseed),1,30,options);
        MC_pars(1,k) = x;
        MC_loss(1,k) = fval;
    end
%    delete(mypool);
    MC_pars_all(smn,:) = MC_pars;
    MC_loss_all(smn,:) = MC_loss;   
end
save MC_pars_all.txt MC_pars_all -ASCII;
save MC_loss_all.txt MC_loss_all -ASCII;
%%
%------Load simulated samples of Gini and Top10, and the actual data
% Note: If "x follows gamma(a,b),then y=1/x follows inverse-gamma(a,1/b)".
sgm_est      = load('optimal_pars.txt');
sgmiv_est    = 1./sgm_est;
MC_sigma_all = load('MC_pars_all.txt');
MC_sgmiv_all = 1./MC_sigma_all;
num_country  = 8;
Nsample      = size(MC_sgmiv_all,1);
PDs          = zeros(num_country,2);
Gtest_sgmiv  = zeros(1,num_country);
SEs_INV      = zeros(num_country,1);  % Standard Error of estimated 1/sigma
SEs_sgm      = zeros(num_country,1);  % Standard Error of estimated sigma
CDF_sgmiv    = zeros(num_country,1);  % CDF value at the estimated 1/sigma
PV_sgmiv     = zeros(num_country,1);  % P-value of the estimated 1/sigma
CTV          = zeros(num_country,1);  % Critical Value for the 95% CI
for i=1:num_country
%------ Test whether the simulated sample of sigma squre follows a gamma    
    [PDs(i,:),CIs9] = gamfit(MC_sgmiv_all(:,i),0.01);
    [~,CIs5]        = gamfit(MC_sgmiv_all(:,i),0.05);
    [~,CIs0]        = gamfit(MC_sgmiv_all(:,i),0.10);
    if (PDs(i,1)>=CIs9(1,1) && PDs(i,1)<=CIs9(2,1)) && (PDs(i,2)>=CIs9(1,2) && PDs(i,2)<=CIs9(2,2))
       Gtest_sgmiv(1,i) = 0.99;
    else
       if (PDs(i,1)>=CIs5(1,1) && PDs(i,1)<=CIs5(2,1)) && (PDs(i,2)>=CIs5(1,2) && PDs(i,2)<=CIs5(2,2)) 
          Gtest_sgmiv(1,i) = 0.95;
       else
          if (PDs(i,1)>=CIs0(1,1) && PDs(i,1)<=CIs0(2,1)) && (PDs(i,2)>=CIs0(1,2) && PDs(i,2)<=CIs0(2,2)) 
             Gtest_sgmiv(1,i) = 0.90;              
          end
       end
    end
%------ Calculate the standard error    
    PDfit  = fitdist(MC_sgmiv_all(:,i),'gamma');
    SD_INV = std(PDfit);
    SEs_INV(i,1) = SD_INV/sqrt(Nsample);
    SD_sgm  = std(MC_sigma_all(:,i));
    SEs_sgm(i,1) = SD_sgm/sqrt(Nsample);
%------ Calculate the P-value for H0: true 1/ = estimated 1/sigma
    CDF_sgmiv(i,1) = gamcdf(sgmiv_est(i,1),PDs(i,1),PDs(i,2));
    PV_sgmiv(i,1) = 1.0-CDF_sgmiv(i,1);
%------ Critical Value for the 95% CI   
    CTV(i,1) = gaminv(0.95,PDs(i,1),PDs(i,2));
end

%------ Save Prb values to the excel file
countries = {'Canada';'China';'France';'Germany';'Italy';'Japan';'UK';'US'};
tableMC2_header = {'Country','Estimated sigma','SE for sigma','P-value for 1/sigma'};
tableMC2 = table(string(countries(:)),sgm_est,SEs_sgm,PV_sgmiv,'VariableNames',matlab.lang.makeValidName(tableMC2_header));
writecell(tableMC2_header,'Generated Tables.xlsx','Sheet','Table4.Estimated sigma','Range','A1:D1');
writetable(tableMC2,'Generated Tables.xlsx','Sheet','Table4.Estimated sigma','Range','A2','WriteVariableNames',false);

%------ Plot distributions of 1/sigma
country_plot1 = 8;                          % Country to be plotted (US)
  gamx1_min = min(MC_sgmiv_all(:,country_plot1));
  gamx1_max = max(MC_sgmiv_all(:,country_plot1));
  gamx1_dif = (gamx1_max-gamx1_min)/(Nsample-1);
  fit1_x = [gamx1_min:gamx1_dif:gamx1_max]';
  fit1_y = gampdf(fit1_x,PDs(country_plot1,1),PDs(country_plot1,2));
  CRT_R1 = gaminv(0.95,PDs(country_plot1,1),PDs(country_plot1,2));
country_plot2 = 2;                          % Country to be plotted (China)
  gamx2_min = min(MC_sgmiv_all(:,country_plot2));
  gamx2_max = max(MC_sgmiv_all(:,country_plot2));
  gamx2_dif = (gamx2_max-gamx2_min)/(Nsample-1);
  fit2_x = [gamx2_min:gamx2_dif:gamx2_max]';
  fit2_y = gampdf(fit2_x,PDs(country_plot2,1),PDs(country_plot2,2));
  CRT_R2 = gaminv(0.95,PDs(country_plot2,1),PDs(country_plot2,2));
ngrid   = 20;
[fp1,xp1]  = ecdf(MC_sgmiv_all(:,country_plot1));
  his1_x   = min(xp1):(max(xp1)-min(xp1))/(ngrid-1):max(xp1);
  his1_y   = ecdfhist(fp1,xp1,ngrid);
  barcolr1 = [.3 .75 .93];
[fp2,xp2]  = ecdf(MC_sgmiv_all(:,country_plot2));
  his2_x   = min(xp2):(max(xp2)-min(xp2))/(ngrid-1):max(xp2);
  his2_y   = ecdfhist(fp2,xp2,ngrid);
  barcolr2 = [1 .4 .4];  
figHandle = create_new_figure('Figure.MC Distribution of sigmas');
h1a=bar(his1_x,his1_y,'FaceColor',barcolr1,'EdgeColor',barcolr1,'LineWidth',1,'FaceAlpha',0.5); hold on;
  h1b=plot(fit1_x,fit1_y,'color','blue','LineWidth',3); hold on;
h2a=bar(his2_x,his2_y,'FaceColor',barcolr2,'EdgeColor',barcolr2,'LineWidth',1,'FaceAlpha',0.5); hold on;
  h2b=plot(fit2_x,fit2_y,'color','red','LineWidth',3); hold on;
dash_yR = gampdf(CRT_R1,PDs(country_plot1,1),PDs(country_plot1,2));
  plot(CRT_R1*ones(101,1),[0:dash_yR/100:dash_yR]','b--','LineWidth',3); hold on; 
  h1c=scatter(sgmiv_est(country_plot1,1),gampdf(sgmiv_est(country_plot1,1),PDs(country_plot1,1),PDs(country_plot1,2)),400,'o','b','MarkerFaceColor','b');
  hold on;
dash_yR = gampdf(CRT_R2,PDs(country_plot2,1),PDs(country_plot2,2));
  plot(CRT_R2*ones(101,1),[0:dash_yR/100:dash_yR]','r--','LineWidth',3); hold on;
  h2c=scatter(sgmiv_est(country_plot2,1),gampdf(sgmiv_est(country_plot2,1),PDs(country_plot2,1),PDs(country_plot2,2)),400,'o','r','MarkerFaceColor','r');
  hold off;
axis_value=axis;
  axis([axis_value(1) axis_value(2) 0 max(max(his1_y),max(his2_y))*1.1]);
  xlabel('Simulated 1/\bf\sigma');    ylabel('Frequency');
text(CRT_R1*0.999,max(fit1_y)*0.25,'CDF=95%','FontName','Times','Fontsize',12,'Color','b','FontWeight','bold');
  text(CRT_R2*0.999,max(fit1_y)*0.25,'CDF=95%','FontName','Times','Fontsize',12,'Color','r','FontWeight','bold');
set(gca,'FontSize',12,'FontName','Times New Roman');
  L=legend([h1a,h1b,h1c,h2a,h2b,h2c],'Histogram of simulated 1/{\bf\sigma} (US)','Fitted Gamma pdf (US)','Estimated 1/{\bf\sigma} (US)','Histogram of simulated 1/{\bf\sigma} (China)','Fitted Gamma pdf (China)','Estimated 1/{\bf\sigma} (China)','location','southoutside','NumColumns',3);
  set(L,'Orientation','horizontal','Box','off','FontSize',12,'FontName','Times New Roman');
set(figHandle,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 12 7]);
%print(figHandle,'-dpng','-r1080','Figure.MC Distribution of sigmas.jpg');


%%        Extension 1: Tax Experiment with Neutral Revenue
%------Replace income tax by wealth tax
%      Step1: Set "Tax_expmt" to 1 on Line 916
%      Step2: Rerun "Section 1" and "Section 3"
%      Step3: Run the following codes
% Table of the average rate in each country
if ~exist('TaxRates','var') || isempty(TaxRates) || ~exist('countries','var') || isempty(countries)
    section3_results = recompute_section3_results_from_settings();
    TaxRates    = section3_results.TaxRates;
    countries   = section3_results.countries;
    wealth_gini = section3_results.wealth_gini;
    t10s        = section3_results.t10s;
    up_rate     = section3_results.up_rate;
    Shorrocks_S = section3_results.Shorrocks_S;
end
Mean_rate = mean(TaxRates,1)';
tableX1a_header = {'Country','Tax Rate'};
tableX1a = table([string(countries(:));"Average"],[Mean_rate;mean(Mean_rate)], ...
    'VariableNames',matlab.lang.makeValidName(tableX1a_header));
writecell(tableX1a_header,'Generated Tables.xlsx','Sheet','Table8a.Wealth Tax','Range','A1:B1');
writetable(tableX1a,'Generated Tables.xlsx','Sheet','Table8a.Wealth Tax','Range','A2','WriteVariableNames',false);
% Table of the estimated Gini and Top share in each country
tableX1b_header = {'Country','Gini','top10','transition','SI'};
tableX1b = table(string(countries(:)),wealth_gini,t10s,up_rate,Shorrocks_S, ...
    'VariableNames',matlab.lang.makeValidName(tableX1b_header));
writecell(tableX1b_header,'Generated Tables.xlsx','Sheet','Table8b.Estimated Gini&Top (W)','Range','A1:E1');
writetable(tableX1b,'Generated Tables.xlsx','Sheet','Table8b.Estimated Gini&Top (W)','Range','A2','WriteVariableNames',false);

%------Replace income tax by consumption tax
%      Step1: Set "Tax_expmt" to 2 on Line 916
%      Step2 and Step3 are same as the above 
%      Note: Don't cover the previous Table8a & Table8b

%%        Extension 2a: Different Assets (US)
% First, run "Section 1" and then directly skip to this section
ext3_country = 8;
  effort = effort_rates(ext3_country,1);  
  growth = growth_rates(ext3_country,1);
  tax = tax_rates(ext3_country,1);
  allowance = allowances(ext3_country,1);
  C0 = consum0(ext3_country,1); 
  C1 = consum1(ext3_country,1);
  act_giniw = act_giniW(ext3_country,1);    
  act_10 = act_top10(ext3_country,1);
  gamma = labor_pwr(ext3_country,1);
  ky_rt = KY_ratio(ext3_country,1);    
options = optimset('Display','none','TolX',1e-3,'TolFun',1e-6);
  [x,fval] = fminbnd(@(x) country_optimal2(x,C0,C1,effort,growth,tax,allowance,act_giniw,act_10,npel,nper,gamma,ky_rt),1,30,options);
  ext3_estimate(1,1) = x;
  ext3_estimate(1,2) = fval;   
save ext3_estimate.txt ext3_estimate -ASCII;

%%         Extension 2a Post-estimation
[equity,~,~] = xlsread('country_data','Diff_Assets_Extension','A2:B6');
sigma_loss   = load('ext3_estimate.txt');
sigma        = sigma_loss(1,1);
npel_aix     = 1:npel;
nper_aix     = (1:nper)';
Lorenz_period= nper;
Y  = npel;       
M1 = 1;
M2 = 1;
M3 = 1;
M4 = 1;
M5 = 1;
wealth_gini = zeros(num_country,1);
t10s      = zeros(num_country,1);
b50s      = zeros(num_country,1);
up_rate   = zeros(num_country,1);
down_rate = zeros(num_country,1);
W_sorted_end = zeros(npel,num_country);
W_sorted_end1= zeros(npel,num_country);
mobility_mat = zeros(Mobl_dim,Mobl_dim,num_country); % Social Mobility matrix (mid to end)
Shorrocks_S  = zeros(num_country,1);
ext3_country = 8;
W0     = ones(npel,1)*KY_ratio(ext3_country,1);
effort = effort_rates(ext3_country,1);
growth = growth_rates(ext3_country,1);
tax    = tax_rates(ext3_country,1);
allowance = allowances(ext3_country,1);
C0 = consum0(ext3_country,1);    
C1 = consum1(ext3_country,1);
gamma = labor_pwr(ext3_country,1);
%-------Call simulation function
gamma_profile = assign_chronic_gamma(gamma, npel, M5);
%-------Call simulation function
[sim_wealth,~,~,RR_other] = fit_sime3(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,effort,sigma,growth,tax,allowance,gamma_profile,equity);
sort_wealth = sort(sim_wealth);   
t10s(ext3_country,1) = sum(sort_wealth(npel*(1-top)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
b50s(ext3_country,1) = sum(sort_wealth(1:npel*TM_btm,nper),1)/sum(sort_wealth(:,nper),1);    
W_sorted_end(:,ext3_country) = sort_wealth(:,nper);
W_sorted_end1(:,ext3_country)= W_sorted_end(:,ext3_country)/mean(W_sorted_end(:,ext3_country)); % Standardized
%-------Lorenz Curve & Gini
scale_sortw  = sort_wealth(:,Lorenz_period);         % Do not Scale    
accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));
wealth_gini(ext3_country,1) = gini_coeff(accum_wealth);     
%-------Calculate upward & downward rates from youth to the old
Threshold_top = sort_wealth(npel*(1-TM_top),nper);   % Lower bound of the top "TM_top"
Threshold_btm = sort_wealth(npel*TM_btm+1,nper);     % Upper bound of the bottom "TM_btm"
sort_mid = sortrows(sim_wealth,mid_per);
sort_mid_nper = sort_mid(:,nper);  % End-period wealth sorted by the mid-period wealth
up_rate(ext3_country,1)   = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
down_rate(ext3_country,1) = sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<Threshold_btm)/(npel*TM_top);    
%-------Social Mobility Matrix ("Mobl_dim" by "Mobl_dim")
bound_btm = zeros(Mobl_dim-1,1);
for i=1:Mobl_dim-1
    bound_btm(i,1) = sort_wealth(npel/Mobl_dim*i,nper);
end
% Mobility matrix [row: bottom to top at the beginning, column: bottom to top at then end]
for j=1:Mobl_dim
    mobility_mat(j,1,ext3_country)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(1,1))/(npel/Mobl_dim);
    for i=2:Mobl_dim-1
        mobility_mat(j,i,ext3_country)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(i-1,1) ...
        & sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(i,1))/(npel/Mobl_dim);
    end
    mobility_mat(j,Mobl_dim,ext3_country)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(Mobl_dim-1,1))/(npel/Mobl_dim); 
end
Shorrocks_S(ext3_country,1) = (Mobl_dim-trace(mobility_mat(:,:,ext3_country)))/(Mobl_dim-1);  % Static Shorrocks Index

%               Estimates of the Thickness Indice
index_start = npel*0.01+1;             % Top range (ie. top1%) started from
index_end   = npel*0.2+1;              % Top range (ie. top20%) ended up
index_every = 1;                       % Agents added in each rolling 
allreg=size((index_start:index_every:index_end),2); % Regressions number
regB  = zeros(2,allreg,num_country);                % Regression coefs
RegSt = zeros(4,allreg,num_country);                % Regression statistics
clear k;
W_sorted_reverse = sort(W_sorted_end1,'descend');   % Reverse sorting
for nreg = 1:allreg
    clear regx regy RegSt0;   
    regy = log((1:index_start-1+(nreg-1)*index_every)');        
    regx(:,1) = ones(index_start-1+(nreg-1)*index_every,1);
    k = 8;
    regx(:,2) = log(W_sorted_reverse(1:index_start-1+(nreg-1)*index_every,k));
    [regB(:,nreg,k),~,~,~,RegSt0] = regress(regy,regx);
    format short g;RegSt(:,nreg,k) = (roundn(RegSt0,-4))';       
end
% Make Tables
Index_X3d(1,1) = -1*regB(2,1,ext3_country);                       % Top 1%
Index_X3d(1,2) = -1*regB(2,round(npel*0.02+1),ext3_country);      % Top 3%
Index_X3d(1,3) = -1*regB(2,round(npel*0.09+1),ext3_country);      % Top 10%
Index_X3d(1,4) = -1*regB(2,round(npel*0.19+1),ext3_country);      % Top 20%
% Save Gini, Top10%, transition Probability, SI, & thickness indices 
table9a_column1={'Gini';'top10';'transition';'SI';'1% Thickness';'3% Thickness';'10% Thickness';'20% Thickness'};
table9a_column2=[wealth_gini(ext3_country,1);t10s(ext3_country,1);up_rate(ext3_country,1);Shorrocks_S(ext3_country,1);Index_X3d'];
table9a_header = {'Metric','Value'};
table9a = table(string(table9a_column1),table9a_column2,'VariableNames',matlab.lang.makeValidName(table9a_header));
writecell(table9a_header,'Generated Tables.xlsx','Sheet','Table9a.Model0_Gini','Range','A1:B1');
writetable(table9a,'Generated Tables.xlsx','Sheet','Table9a.Model0_Gini','Range','A2','WriteVariableNames',false);
% Save implied return rate of non-financial assets
table9b_header = {'Description','Value'};
rr_other_values = RR_other(2:end,1);
additional_rows = max(numel(rr_other_values)-1,0);
table9b_descriptions = [{'Average'};{'Return rate of other assets'};repmat({''},additional_rows,1)];
table9b_values = [mean(rr_other_values);rr_other_values];
table9b = table(string(table9b_descriptions),table9b_values,'VariableNames',matlab.lang.makeValidName(table9b_header));
writecell(table9b_header,'Generated Tables.xlsx','Sheet','Table9b.Extension3_RR_other','Range','A1:B1');
writetable(table9b,'Generated Tables.xlsx','Sheet','Table9b.Extension3_RR_other','Range','A2','WriteVariableNames',false);


%%         Extension 2b: Model1— Raise return rate of equity
% Step1: Multiplies the equity return rate in the data file by 1+10%
% Step2: Rerun "Extension 2a Post-estimation" above
% Note: Don't cover the previous results

%%         Extension 2c: Model2— Raise equity ratio of the poor
% Step1: Set equity asset share of each group equal to that of the rich
% Step2: Rerun "Extension 2a Post-estimation" above
% Note: Don't cover the previous results


%%                             Functions
%==========================================================================
% Function 1: Estimate parameters across countries                         
function loss = country_optimal(sigma,C0,C1,beta1,growth,tax,allowance,actginiw,act10,npel,nper,gamma,ky_rt,tmseed)
% objetive loss function to estimate the parameters
W0 = ones(npel,1)*ky_rt;   % Initial wealth endowment
M1 = 1;  % Effort
M2 = 1;  % Growth
M3 = 1;  % Tax
M4 = 1;  % Capital
M5 = 1;  % Rent privilege scenario (disabled during estimation)
Y = npel;  % Output
gamma_profile = assign_chronic_gamma(gamma, npel, M5);
%-------Call simulation function
[sim_wealth,~,~,~] = fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance,gamma_profile,tmseed);
sort_wealth = sort(sim_wealth); 
Lorenz_period = nper; % Set the period when we plot Lorenz & calculate gini   
%-------Calculate Inequality  
scale_sortw = sort_wealth(:,Lorenz_period); % Do not Scale
accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));    
gini_w = gini_coeff(accum_wealth);        
top10 = sum(sort_wealth(npel*(1-0.1)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
%-------Calculate Loss
loss = (gini_w-actginiw)^2+(top10-act10)^2;    
end

%==========================================================================
% Function 2: Simulate data using the Single-factor Model/the Hybrid Model                         
function[WE,ICM,ICM2,TaxRate]=fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance,gamma,tmseed)
% NEW INPUT----------------------------------------------------------------
% M1:Labour model, M2:Growth model, M3:Tax model, M4:Capital model, M5:Rent
%     privileges indicator
% Mi can take on two values, "1" means being considered while "0" means not

% W0:        Initial wealth at the begining of period 1
% nper:      Number of periods
% c0:        Individual constant consumption level (lower bound)
% c1:        Proportional consumption out of individual income
% edr:       Mrginal return in terms of human capital
% growth:    A net exogenous growth rate
% cpt:       Consumption tax rate
% kpt:       Capital tax rate
% OUTPUT ------------------------------------------------------------------
% WE:        Individual Wealth at the end of each period (npel by nper)
% ICM:       Individual total income in each period
npel = size(W0,1);         % npel: Number of people
WB   = zeros(npel,nper+1); % WB: Wealth at the beginning of each period
WE   = zeros(npel,nper);
CSM  = zeros(npel,nper); % Individual consumption if consumption tax is in
CSM0 = zeros(npel,nper); % Individual consumption if consumption tax is out
ICM  = zeros(npel,nper); % Individual income before any subsidy or transfer
ICM2 = zeros(npel,nper); % Individual total income including transfer
WB(:,1)   = W0(:,1);
rawprob   = zeros(npel,1);
inv_share = zeros(npel,1);

M5 = logical(M5);
rent_mask = false(npel,1);
if M5
    num_rent_seekers = max(1, ceil(0.01 * npel));
    rent_mask(1:num_rent_seekers,1) = true;
else
    num_rent_seekers = 0;
end
TaxRate   = zeros(nper,1);

if isscalar(gamma)
    labor_weights = repmat(gamma, npel, 1);
elseif numel(gamma) == npel
    labor_weights = reshape(gamma, [], 1);
else
    error('Code_Hybrid_Model:InvalidGammaProfile', ...
        'Gamma must be a scalar or match the number of agents.');
end

labor_weights = max(labor_weights, 0);
sum_labor_weights = sum(labor_weights);
avg_labor_share = sum_labor_weights / npel;
capital_weights = max(1 - labor_weights, 0);
sum_capital_weights = sum(capital_weights);

Tax_expmt = 0;  % Taxation experiment trigger,
          % 1: income tax is replaced by wealth tax; 2: by consumption tax

rng(tmseed);  % Fix the time seed for random generator (this syntax is newest)
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
    if M3==1
       T = tax;
    else
       T = 0.0;
    end
    period_output = Y*(1+G)^(j-1);
    labor_total = period_output * avg_labor_share;
    if sum_labor_weights>0
        labor_unit = labor_total / sum_labor_weights;
    else
        labor_unit = 0;
    end
    capital_share_factor = max(1 - avg_labor_share, 0);
    capital_total = period_output * capital_share_factor;
    if sum_capital_weights>0
        capital_unit = capital_total / sum_capital_weights;
    else
        capital_unit = 0;
    end
% Creates probabilities determined by Human Capital
    rawprob(:,1) = 1 - 1./(ones(npel,1)+exp(ones(npel,1)*log(1/(npel-1))+beta1*effort));
    prob1 = rawprob/sum(rawprob); % Logit probability
% Creates probabilities determined by Physical Capital Investment
    if j==1
       inv_share(:,1) = 1/npel;          % No investment effect in period 1
    else
       inv_share(:,1) = WE(:,j-1)/sum(WE(:,j-1),1);
    end
%------2nd, Distribute GDP sequentially
    y1 = zeros(npel,1);                       % Initial income distribution
    if M4==0
       if labor_unit>0 && sum_labor_weights>0
          gain0 = randsrc(npel,1,[1:1:npel; prob1']);       % Probability of labour income
          for i=1:npel
              idx = gain0(i,1);
              y1(idx,1) = y1(idx,1) + labor_unit*labor_weights(idx,1);
          end
       end
    elseif M4==1 && M1==0
       if capital_unit>0 && sum_capital_weights>0
          gain1 = randsrc(npel,1,[1:1:npel; inv_share']); % Probability of capital income
          for i=1:npel
              idx = gain1(i,1);
              y1(idx,1) = y1(idx,1) + capital_unit*capital_weights(idx,1);
          end
       end
    elseif M1==1 && M4==1
       if labor_unit>0 && sum_labor_weights>0
          gain0 = randsrc(npel,1,[1:1:npel; prob1']);     % Probability of labour income
          for i=1:npel                        % 1st, distribute labour income
              idx = gain0(i,1);
              y1(idx,1) = y1(idx,1) + labor_unit*labor_weights(idx,1);
          end
       end
       if capital_unit>0 && sum_capital_weights>0
          gain1 = randsrc(npel,1,[1:1:npel; inv_share']); % Probability of capital income
          for i=1:npel                        % 2nd, distribute capital income
              idx = gain1(i,1);
              y1(idx,1) = y1(idx,1) + capital_unit*capital_weights(idx,1);
          end
       end
    end
%------3rd, Calculate income tax burden rate if there is no other taxes
    taxable_y = zeros(npel,1);
    for i=1:npel
        if y1(i,1)>allowance*(Y*(1+G)^(j-1)/npel)
           taxable_y(i,1) = y1(i,1) - allowance*(Y*(1+G)^(j-1)/npel);                
        end
    end
    clear taxpay_y0;
    taxpay_y0 = taxable_y*T;
    if any(rent_mask)
        taxpay_y0(rent_mask,1) = 0;
    end
    T_burden = sum(taxpay_y0);
%=======Case1: Only income tax    
    if Tax_expmt==0
       taxpay_y = taxpay_y0;
       Taxall_y = T_burden;
     % Subsidy to the very poor
       y2 = zeros(npel,1);          % 1st time transfer
       for i=1:npel
           if WB(i,j) + y1(i,1) - taxpay_y(i,1) < C0*(1+G)^(j-1)
              y2(i,1) = C0*(1+G)^(j-1);
           end
       end
       if any(rent_mask)
           y2(rent_mask,1) = C0*(1+G)^(j-1);
       end
       num_sub1 = sum(y2(:,1));     % Total amount of the 1st time transfer
       ICM2(:,j) = y1 - taxpay_y;   % Income before transfer
     % Transfer the rest income tax revenue to agents identically
       y3 = zeros(npel,1);          % 2nd time transfer
       for i=1:npel
           if Taxall_y>num_sub1
              y3(i,1) = (Taxall_y-num_sub1)/npel;              
           end
       end
       ICM(:,j) = ICM2(:,j) + y2 + y3;               % Income after 2nd transfer       
     % Consumption according to "ICM2"
       for i=1:npel
           if WB(i,j) + ICM2(i,j) < C0*(1+G)^(j-1) + C1*ICM2(i,j)
              CSM(i,j) = WB(i,j) + ICM2(i,j);
           else
              CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
           end
       end
       if any(rent_mask)
          CSM(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
       end
     % Calculate individual wealth after consumption
       WE(:,j) = WB(:,j) + ICM(:,j) - CSM(:,j);
       TaxRate(j,1) = T;       
%=======Case2: Only wealth tax but with same tax burden as income tax 
    elseif Tax_expmt==1       
       TW_multp = 1;       
       WB_b1 = WB(:,j) + y1;
       Wmean = mean(WB_b1);
      % Wmedian = median(WB_b1);
       taxable_k = zeros(npel,1);           % Taxable wealth for each agent
       for i=1:npel
           if WB_b1(i,1) > TW_multp*Wmean
              taxable_k(i,1) = WB_b1(i,1) - TW_multp*Wmean;
          % if WB_b1(i,1)>TW_multp*Wmedian
          %    taxable_k(i,1) = WB_b1(i,1)-TW_multp*Wmedian;
          % if WB_b1(i,1)>TW_multp*Wtop(round((1.0-0.05)*npel),1)
          %    taxable_k(i,1) = WB_b1(i,1)-TW_multp*Wtop(round((1.0-0.05)*npel),1);
           end
       end
       if any(rent_mask)
           taxable_k(rent_mask,1) = 0;
       end
       denom_taxable = sum(taxable_k(:,1));
       if denom_taxable == 0
           tau_k = 0;
       else
           tau_k = T_burden/denom_taxable;
       end
       TaxRate(j,1) = tau_k;
       clear taxpay_k;
       taxpay_k = taxable_k*tau_k;
       if any(rent_mask)
           taxpay_k(rent_mask,1) = 0;
       end
       Taxall_k = sum(taxpay_k(:,1));
       WB_b2 = WB_b1 - taxpay_k;
     % Subsidy to the very poor
       y2 = zeros(npel,1);
       for i=1:npel
           if WB_b2(i,1)<C0*(1+G)^(j-1)
              y2(i,1) = C0*(1+G)^(j-1);
           end
       end
       if any(rent_mask)
          y2(rent_mask,1) = C0*(1+G)^(j-1);
       end
       num_sub2 = sum(y2);               % Total amount of the 1st transfer
       WB_b3 = WB_b2 + y2;               % Wealth after 1st transfer
       ICM2(:,j) = y1;                   % Income before transfer       
     % Transfer the rest income tax revenue to agents identically
       y3 = zeros(npel,1);               % 2nd time transfer
       for i=1:npel
           if Taxall_k>num_sub2
              y3(i,1) = (Taxall_k-num_sub2)/npel;
           end
       end
       WB_b4 = WB_b3 + y3;                      % Wealth after 2nd transfer        
       ICM(:,j) = ICM2(:,j) + y2 + y3;          % Income after 2nd transfer       
     % Consumption according to "ICM2" and "WB_b3"
       for i=1:npel
           if WB_b3(i,1)<C0*(1+G)^(j-1)+C1*ICM2(i,j)
              CSM(i,j) = WB_b3(i,1);
           else
              CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
           end
       end
       if any(rent_mask)
          CSM(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
       end
     % Calculate individual wealth after consumption
       WE(:,j) = WB_b4 - CSM(:,j);
%=======Case3: Only consumption tax but with same tax burden as income tax 
    elseif Tax_expmt==2             
     % Note: We have to transfer subsidy by public deficit first and then
     %       balance the deficit by consumption tax revenue in thise case.
       WB_c1 = WB(:,j) + y1;
     % Subsidy to the very poor
       y2 = zeros(npel,1);
       for i=1:npel
           if WB_c1(i,1)<C0*(1+G)^(j-1)
              y2(i,1) = C0*(1+G)^(j-1);
           end
       end
       if any(rent_mask)
          y2(rent_mask,1) = C0*(1+G)^(j-1);
       end
       num_sub3 = sum(y2);               % Total amount of the 1st transfer        
       WB_c2 = WB_c1 + y2;               % Wealth after 1st transfer
       ICM2(:,j) = y1;                   % Income after before transfer       
%----- Following steps aim to compute the proper consumption tax rate        
       % Step1. Calculate initial tau_c given the income tax burden
       %       (Assume agents consume following the consumption equation)
       tau_c0 = T_burden/sum(C1*ICM2(:,j),1);
       Tburden_tol = 0.01;  % Tolerance level of solving for tau_c
       % Step2. Calculate the real consumption tax burden given tau_c0
       for i=1:npel
           if WB_c2(i,1) < C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c0)
              CSM(i,j) = WB_c2(i,1);           % Consumption with tau_c0
           else
              CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c0);
           end
           if WB_c2(i,1) < C0*(1+G)^(j-1) + C1*ICM2(i,j)
              CSM0(i,j) = WB_c2(i,1);          % Consumption without tau_c0
           else
              CSM0(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
           end
       end
       if any(rent_mask)
          CSM0(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
          CSM(rent_mask,j) = CSM0(rent_mask,j);
       end
       Tburden_diff = (sum(CSM(:,j),1)-sum(CSM0(:,j),1))/T_burden;
       tau_c = tau_c0;
       while abs(Tburden_diff-1.0) > Tburden_tol
           if Tburden_diff-1.0>Tburden_tol             % tau_c is too large
              tau_c = (tau_c0+tau_c)/2;
           elseif Tburden_diff-1.0<-Tburden_tol        % tau_c is too small
              tau_c = tau_c*(1+0.1);
           end
           for i=1:npel      % Recompute consumption tax burden given tau_c
               if WB_c2(i,1) < C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c)
                  CSM(i,j) = WB_c2(i,1);
               else
                  CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c);
               end
               if WB_c2(i,1) < C0*(1+G)^(j-1)+ C1*ICM2(i,j)
                  CSM0(i,j) = WB_c2(i,1);
               else
                  CSM0(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
               end
           end
           if any(rent_mask)
              CSM0(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
              CSM(rent_mask,j) = CSM0(rent_mask,j);
           end
           Tburden_diff = (sum(CSM(:,j),1)-sum(CSM0(:,j),1))/T_burden;
       end
       % Now the real tau_c has been found!
       TaxRate(j,1) = tau_c;
       for i=1:npel    % Recompute real consumption tax burden given a new tau_c
           if WB_c2(i,1) < C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c)
              CSM(i,j) = WB_c2(i,1);
           else
              CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j)*(1+tau_c);
           end
           if WB_c2(i,1) < C0*(1+G)^(j-1) + C1*ICM2(i,j)
              CSM0(i,j) = WB_c2(i,1);
           else
              CSM0(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
           end
       end
       if any(rent_mask)
          CSM0(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
          CSM(rent_mask,j) = CSM0(rent_mask,j);
       end
     % Transfer the rest consumption revenue to agents identically 
       Taxall_c = sum(CSM(:,j),1)-sum(CSM0(:,j),1);     % Total cpt revenue        
       y3 = zeros(npel,1);                              % 2nd time transfer
       if Taxall_c>num_sub3
          for i=1:npel
              y3(i,1) = (Taxall_c-num_sub3)/npel;
          end
       end
       ICM(:,j) = ICM2(:,j) + y2 + y3;          % Income after 2nd transfer       
       WE(:,j) = WB(:,j) + ICM(:,j) - CSM(:,j);         
    end
    WB(:,j+1) = WE(:,j);
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

function gamma_profile = assign_chronic_gamma(gamma, npel, M5)
%ASSIGN_CHRONIC_GAMMA Apply the chronic gamma adjustment to the agent vector.
%   GAMMA_PROFILE = ASSIGN_CHRONIC_GAMMA(GAMMA, NPEL, M5) expands the scalar
%   GAMMA to an NPEL-by-1 vector. When M5 equals 1, the first one percent of
%   agents are flagged as rent-seekers and receive a gamma that is ten times
%   the average gamma prior to the adjustment. If GAMMA is already a vector
%   with length NPEL it is returned as a column vector before the rent-seeker
%   modification is applied.

if nargin < 3
    M5 = 0;
end

if isscalar(gamma)
    gamma_profile = repmat(gamma, npel, 1);
elseif numel(gamma) == npel
    gamma_profile = reshape(gamma, [], 1);
else
    error('Code_Hybrid_Model:InvalidGammaProfile', ...
        'Gamma must be a scalar or a vector with one entry per agent.');
end

if M5 == 1
    num_rent_seekers = max(1, ceil(0.01 * npel));
    base_mean_gamma = mean(gamma_profile);
    gamma_profile(1:num_rent_seekers,1) = base_mean_gamma * 10;
end

end

%==========================================================================
% Function 1b: Estimate Extension 3 with Different Assets
function loss = country_optimal2(sigma,C0,C1,beta1,growth,tax,allowance,actginiw,act10,npel,nper,gamma,ky_rt)
% objetive loss function to estimate the parameters
W0 = ones(npel,1)*ky_rt;  % Initial wealth endowment
[equity,~,~] = xlsread('country_data','Diff_Assets_Extension','A2:B6');
% equity(1,1) is 'Equity Ratio of Assets for Top1% wealth holder'
% equity(1,2) is 'Equity Ratio of Assets for 90%-99% wealth holder'
% equity(1,3) is 'Equity Ratio of Assets for 50%-90% wealth holder'
% equity(1,4) is 'Equity Ratio of Assets for Bottom50% wealth holder'
% equity(1,5) is  'Equity Return Rate'
M1 = 1;  % Effort 
M2 = 1;  % Growth
M3 = 1;  % Tax
M4 = 1;  % Capital
M5 = 1;  % Rent privilege scenario (disabled during estimation)
Y = npel;  % Output
%-------Call simulation function
[sim_wealth,~,~,~] = fit_sime3(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance,assign_chronic_gamma(gamma,npel,M5),equity);
sort_wealth = sort(sim_wealth); 
Lorenz_period = nper; % Set the period when we plot Lorenz & calculate gini   
%-------Calculate Inequality  
scale_sortw = sort_wealth(:,Lorenz_period); % Do not Scale
accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));    
gini_w = gini_coeff(accum_wealth);        
top10 = sum(sort_wealth(npel*(1-0.1)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
%-------Calculate Loss
loss = (gini_w-actginiw)^2+(top10-act10)^2;    
end

%==========================================================================
% Function 2b: Simulate data using the Single-factor Model/the Hybrid Model                         
function[WE,ICM,ICM2,RR_other]=fit_sime3(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,beta1,sigma,growth,tax,allowance,gamma,equity)
% INPUT----------------------------------------------------------------
% M1:Labour model, M2:Growth model, M3:Tax model, M4:Capital model, M5:Rent
%     privileges indicator
% Mi can take on two values, "1" means being considered while "0" means not
% W0:        Initial wealth at the begining of period 1
% nper:      Number of periods
% c0:        Individual constant consumption level (lower bound)
% c1:        Proportional consumption out of individual income
% edr:       Mrginal return in terms of human capital
% growth:    A net exogenous growth rate
% cpt:       Consumption tax rate
% kpt:       Capital tax rate
% OUTPUT ------------------------------------------------------------------
% WE:        Individual Wealth at the end of each period (npel by nper)
% ICM:       Individual total income in each period
npel = size(W0,1);
WB   = zeros(npel,nper+1);
WE   = zeros(npel,nper);
CSM  = zeros(npel,nper);
ICM  = zeros(npel,nper);
ICM2 = zeros(npel,nper);
WB(:,1)   = W0(:,1);
rawprob   = zeros(npel,1);
inv_share = zeros(npel,1);

M5 = logical(M5);
rent_mask = false(npel,1);
if M5
    num_rent_seekers = max(1, ceil(0.01 * npel));
    rent_mask(1:num_rent_seekers,1) = true;
else
    num_rent_seekers = 0;
end

if isscalar(gamma)
    labor_weights = repmat(gamma, npel, 1);
elseif numel(gamma) == npel
    labor_weights = reshape(gamma, [], 1);
else
    error('Code_Hybrid_Model:InvalidGammaProfile', ...
        'Gamma must be a scalar or match the number of agents.');
end

labor_weights = max(labor_weights, 0);
sum_labor_weights = sum(labor_weights);
avg_labor_share = sum_labor_weights / npel;
capital_weights = max(1 - labor_weights, 0);
sum_capital_weights = sum(capital_weights);

RR_other  = zeros(nper+1,1);     % Return rate of non-financial assets
RR_equity = equity(5,1);         % Equity's Return Rate
WB_equity = zeros(npel,nper+1);  % Individual Equity at the beginning
WB_other  = zeros(npel,nper+1);  % Individual Other Assets at the beginning

rng(1);    % Fix the time seed for random generator (this syntax is newest)
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
equity_share1  = 0.15;                    % Financial Assets at period 1 (Group average)
WB_equity(:,1) = WB(:,1)*equity_share1;   % Financial Assets
WB_other(:,1)  = WB(:,1)-WB_equity(:,1);  % Non-financial Assets9
RR_other(1,1)  = ((1-avg_labor_share)*Y-RR_equity*sum(WB_equity(:,1)))/sum(WB_other(:,1));
for j=1:nper
    if M3==1
       T = tax;
    else
       T = 0;
    end
% ======== New Setting for Extension 3
    y1 = zeros(npel,1);
    % Return of Non-financial assets
    period_output = Y*(1+G)^(j-1);
    labor_total = period_output * avg_labor_share;
    if sum_labor_weights>0
        labor_unit = labor_total / sum_labor_weights;
    else
        labor_unit = 0;
    end
    capital_share_factor = max(1 - avg_labor_share, 0);
    capital_total = period_output * capital_share_factor;
    if sum_capital_weights>0
        capital_unit = capital_total / sum_capital_weights;
    else
        capital_unit = 0;
    end
% Creates probabilities determined by Human Capital
    rawprob(:,1) = 1 - 1./(ones(npel,1)+exp(ones(npel,1)*log(1/(npel-1))+beta1*effort));
    prob1 = rawprob/sum(rawprob);        % Logit probability
% Creates probabilities determined by Physical Capital Investment
    if j==1
       inv_share(:,1) = 1/npel;          % No investment effect in period 1
    else   
       inv_share(:,1) = WB_equity(:,j)/sum(WB_equity(:,j),1);
    end
%------2nd, Distribute the GDP sequentially
    if M4==0
       if labor_unit>0 && sum_labor_weights>0
          gain0 = randsrc(npel,1,[1:npel; prob1']); % Probability of labor income
          for i=1:npel
              idx = gain0(i,1);
              y1(idx,1) = y1(idx,1) + labor_unit*labor_weights(idx,1);
          end
       end
    elseif M4==1 && M1==0
       if capital_unit>0 && sum_capital_weights>0
          gain1 = randsrc(npel,1,[1:npel; inv_share']); % Probability of capital income
          for i=1:npel
              idx = gain1(i,1);
              y1(idx,1) = y1(idx,1) + capital_unit*capital_weights(idx,1);
          end
       end
    elseif M1==1 && M4==1
       if labor_unit>0 && sum_labor_weights>0
          gain0    = randsrc(npel,1,[1:npel; prob1']);  % Probability of labour income
          for i=1:npel                         % 1st, distribute labour income
              idx = gain0(i,1);
              y1(idx,1) = y1(idx,1) + labor_unit*labor_weights(idx,1);
          end
       end
       y1_NFA   = WB_other(:,j)*RR_other(j,1);
       y1a   = y1 + y1_NFA;                 % 2nd, distribute other returns
       gain1 = randsrc(npel,1,[1:npel; inv_share']); % Probability of equity return
       y1b   = zeros(npel,1);               % 3rd, distribute equity return
       WB_equity_dist = sum(WB_equity(:,j)*RR_equity,1);
       for i=1:npel
           y1b(gain1(i,1),1) = y1b(gain1(i,1),1) + WB_equity_dist/npel;
       end
       y1 = y1a + y1b;
    end
% ======== New Setting Ends
%------ Charge tax
    taxpay = zeros(npel,1);
    for i=1:npel
        if y1(i,1)>allowance*(Y*(1+G)^(j-1)/npel)
           taxpay(i,1) = (y1(i,1)-allowance*(Y*(1+G)^(j-1)/npel))*T;
        end
    end
    if any(rent_mask)
        taxpay(rent_mask,1) = 0;
    end
    Taxall = sum(taxpay);
%------ Subsidy to the very poor
    y2 = zeros(npel,1);
    for i=1:npel
        if WB(i,j)+ y1(i,1)- taxpay(i,1) < C0*(1+G)^(j-1)
           y2(i,1) = C0*(1+G)^(j-1);
        end
    end
    if any(rent_mask)
        y2(rent_mask,1) = C0*(1+G)^(j-1);
    end
       num_sub = sum(y2);           % Total amount of the 1st time transfer
       ICM2(:,j) = y1 - taxpay;     % Income before transfer  
     % Transfer the rest income tax revenue to agents identically
       y3 = zeros(npel,1);          % 2nd time transfer
       for i=1:npel
           if Taxall>num_sub
              y3(i,1) = (Taxall-num_sub)/npel;              
           end
       end
       ICM(:,j) = ICM2(:,j) + y2 + y3;          % Income after 2nd transfer  
     % Consumption according to "ICM2"
       for i=1:npel
           if WB(i,j) + ICM2(i,j) < C0*(1+G)^(j-1) + C1*ICM2(i,j)
              CSM(i,j) = WB(i,j) + ICM2(i,j);
           else
              CSM(i,j) = C0*(1+G)^(j-1) + C1*ICM2(i,j);
           end
       end
       if any(rent_mask)
          CSM(rent_mask,j) = C0*(1+G)^(j-1) + C1*ICM2(rent_mask,j);
       end
%   Calculate individual wealth after consumption
    WE(:,j)  = WB(:,j) + ICM(:,j) - CSM(:,j);
    WB(:,j+1)= WE(:,j);
%------Introduce 2 Types of Assets and Different Return Rates
%------Step1.Confirm individual asset components
    Wsort_temp = sort(WE(:,j),1,'descend');
    Equity_ratio = zeros(npel,1);
    for i=1:npel
        if WE(i,j)>=Wsort_temp(round(npel*0.01),1)
           Equity_ratio(i,1) = equity(1,1);
        elseif WE(i,j)>=Wsort_temp(round(npel*0.1),1) && WE(i,j)<Wsort_temp(round(npel*0.01),1)
           Equity_ratio(i,1) = equity(2,1);
        elseif WE(i,j)>=Wsort_temp(round(npel*0.5),1) && WE(i,j)<Wsort_temp(round(npel*0.1),1)
           Equity_ratio(i,1) = equity(3,1);
        else
           Equity_ratio(i,1) = equity(4,1);           
        end
    end
    WB_equity(:,j+1) = WE(:,j).*Equity_ratio;        % Financial Assets
    WB_other(:,j+1)  = WE(:,j)-WB_equity(:,j+1);     % Non-financial Assets
    WB_equity_sum = sum(WB_equity(:,j+1),1);
    WB_other_sum  = sum(WB_other(:,j+1),1);
%------Step2.Calculate the return rate of other assets
    % Since (1-alpha)Y_t+1 = RR_equity*equity+RR_other*other_assets, we have
RR_other(j+1,1) = ((1-avg_labor_share)*(Y*(1+G)^j)-RR_equity*WB_equity_sum)/WB_other_sum;
end

end


function results = recompute_section3_results_from_settings()
settings_file = 'settings_hybrid.mat';
if exist(settings_file,'file') ~= 2
    error('Code_Hybrid_Model:MissingSettings', ...
        ['Required simulation settings are unavailable. ', ...
        'Please run Section 1 before executing Extension 1.']);
end
settings_data = load(settings_file);
required_fields = {'npel','nper','mid_per','num_country','top','TM_top','TM_btm','Mobl_dim'};
for iField = 1:numel(required_fields)
    if ~isfield(settings_data, required_fields{iField})
        error('Code_Hybrid_Model:IncompleteSettings', ...
            ['The settings file is missing the "', required_fields{iField}, '" parameter. ', ...
            'Please rerun Section 1 to regenerate the settings.']);
    end
end
npel        = settings_data.npel;
nper        = settings_data.nper;
mid_per     = settings_data.mid_per;
num_country = settings_data.num_country;
top         = settings_data.top;
TM_top      = settings_data.TM_top;
TM_btm      = settings_data.TM_btm;
Mobl_dim    = settings_data.Mobl_dim;

[actual_data,txt,~] = xlsread('country_data','load_data','A1:L9');
if size(actual_data,1) < num_country
    error('Code_Hybrid_Model:InsufficientCountryData', ...
        'The input data file does not contain enough country entries.');
end
sigmas       = load('optimal_pars.txt');
effort_rates = actual_data(:,1);
growth_rates = actual_data(:,2);
allowances   = actual_data(:,3);
tax_rates    = actual_data(:,4);
consum0      = actual_data(:,5);
consum1      = actual_data(:,6);
labor_pwr    = actual_data(:,9);
KY_ratio     = 3.0*ones(num_country,1);
tmseed       = 1;
if isfield(settings_data,'M5')
    M5 = settings_data.M5;
else
    M5 = 0;
end

[TaxRates, W_sorted_end, W_sorted_end1, mobility_mat, wealth_gini, t10s, ...
    b50s, up_rate, down_rate, Shorrocks_S] = compute_section3_simulations_from_data( ...
    npel, nper, mid_per, num_country, top, TM_top, TM_btm, Mobl_dim, ...
    KY_ratio, effort_rates, growth_rates, allowances, tax_rates, ...
    consum0, consum1, sigmas, labor_pwr, tmseed, M5);

results = struct('TaxRates', TaxRates, ...
    'countries', {txt(2:1+num_country,1)}, ...
    'wealth_gini', wealth_gini, ...
    't10s', t10s, ...
    'up_rate', up_rate, ...
    'Shorrocks_S', Shorrocks_S, ...
    'b50s', b50s, ...
    'down_rate', down_rate, ...
    'mobility_mat', mobility_mat, ...
    'W_sorted_end', W_sorted_end, ...
    'W_sorted_end1', W_sorted_end1);
end


function data = read_sheet_with_fallback(workbook, sheet, range)
%READ_SHEET_WITH_FALLBACK Read worksheet data using multiple strategies.
%   DATA = READ_SHEET_WITH_FALLBACK(WORKBOOK, SHEET, RANGE) attempts to read
%   the specified RANGE from SHEET in WORKBOOK using a sequence of
%   increasingly broad import strategies. This helps accommodate MATLAB
%   environments where 'basic' mode is unavailable or sheet name matching is
%   sensitive to trailing spaces.

last_error = [];

% First try the legacy xlsread basic mode, which avoids requiring Excel on
% Windows but is not available in all environments.
try
    data = xlsread(workbook, sheet, range, 'basic');
    if ~isempty(data)
        return
    end
catch err
    last_error = err;
end

% Fall back to standard xlsread behaviour.
try
    data = xlsread(workbook, sheet, range);
    if ~isempty(data)
        return
    end
catch err
    last_error = err;
end

% Use readmatrix when available (introduced in R2019a) as an additional
% compatibility option.
if exist('readmatrix', 'file') == 2
    try
        data = readmatrix(workbook, 'Sheet', sheet, 'Range', range);
        return
    catch err
        last_error = err;
    end
end

% If all strategies fail, propagate the most recent error.
if ~isempty(last_error)
    rethrow(last_error);
end

% As a final safeguard, throw a generic error if we never captured a
% specific failure reason (which should not normally occur).
if isnumeric(sheet)
    sheet_label = num2str(sheet);
else
    sheet_label = sheet;
end
error('Code_Hybrid_Model:WorksheetReadFailed', ...
    'Unable to read worksheet "%s" from %s.', sheet_label, workbook);
end


function normalized = normalize_sheet_name(name)
%NORMALIZE_SHEET_NAME Generate a canonical representation of a sheet name.
%   NORMALIZED = NORMALIZE_SHEET_NAME(NAME) lower-cases NAME, removes leading
%   and trailing whitespace, and strips punctuation characters that commonly
%   vary across spreadsheet exports. This helps us perform resilient
%   comparisons when deciding which worksheet to load.

if isstring(name)
    name = char(name);
end

if ~ischar(name)
    error('Code_Hybrid_Model:InvalidSheetName', ...
        'Sheet names must be provided as character vectors or strings.');
end

normalized = lower(strtrim(name));
normalized = regexprep(normalized, '[\s\._-]+', ' ');
normalized = regexprep(normalized, '\s+', '');
normalized = regexprep(normalized, '[^a-z0-9]', '');
end


function [TaxRates, W_sorted_end, W_sorted_end1, mobility_mat, wealth_gini, t10s, ...
    b50s, up_rate, down_rate, Shorrocks_S] = compute_section3_simulations_from_data( ...
    npel, nper, mid_per, num_country, top, TM_top, TM_btm, Mobl_dim, KY_ratio, ...
    effort_rates, growth_rates, allowances, tax_rates, consum0, consum1, ...
    sigmas, labor_pwr, tmseed, M5)

Y  = npel;
M1 = 1;
M2 = 1;
M3 = 1;
M4 = 1;
if nargin < 19, M5 = 0; end;
M5 = logical(M5);
Lorenz_period= nper;
wealth_gini  = zeros(num_country,1);
t10s         = zeros(num_country,1);
b50s         = zeros(num_country,1);
up_rate      = zeros(num_country,1);
down_rate    = zeros(num_country,1);
W_sorted_end = zeros(npel,num_country);
W_sorted_end1= zeros(npel,num_country);
mobility_mat = zeros(Mobl_dim,Mobl_dim,num_country);
Shorrocks_S  = zeros(num_country,1);
TaxRates     = zeros(nper,num_country);
for k=1:num_country
    W0 = ones(npel,1)*KY_ratio(k,1);
    effort = effort_rates(k,1);
    growth = growth_rates(k,1);
    tax = tax_rates(k,1);
    allowance = allowances(k,1);
    C0 = consum0(k,1);
    C1 = consum1(k,1);
    sigma = sigmas(k,1);
    gamma = labor_pwr(k,1);
    gamma_profile = assign_chronic_gamma(gamma, npel, M5);
    [sim_wealth,~,~,TaxRate] = fit_sim(W0,M1,M2,M3,M4,M5,Y,nper,C0,C1,effort,sigma,growth,tax,allowance,gamma_profile,tmseed);
    sort_wealth   = sort(sim_wealth);
    TaxRates(:,k) = TaxRate;
    t10s(k,1) = sum(sort_wealth(npel*(1-top)+1:npel,nper),1)/sum(sort_wealth(:,nper),1);
    b50s(k,1) = sum(sort_wealth(1:npel*TM_btm,nper),1)/sum(sort_wealth(:,nper),1);
    W_sorted_end(:,k) = sort_wealth(:,nper);
    W_sorted_end1(:,k)= W_sorted_end(:,k)/mean(W_sorted_end(:,k));
    scale_sortw  = sort_wealth(:,Lorenz_period);
    accum_wealth = cumsum(scale_sortw(:,1))/sum(scale_sortw(:,1));
    wealth_gini(k,1) = gini_coeff(accum_wealth);
    Threshold_top = sort_wealth(npel*(1-TM_top),nper);
    Threshold_btm = sort_wealth(npel*TM_btm+1,nper);
    sort_mid      = sortrows(sim_wealth,mid_per);
    sort_mid_nper = sort_mid(:,nper);
    up_rate(k,1)  = sum(sort_mid_nper(1:npel*TM_btm,1)>Threshold_top)/(npel*TM_btm);
    down_rate(k,1)= sum(sort_mid_nper(npel*(1-TM_top)+1:npel,1)<Threshold_btm)/(npel*TM_top);
    bound_btm = zeros(Mobl_dim-1,1);
    for i=1:Mobl_dim-1
        bound_btm(i,1) = sort_wealth(npel/Mobl_dim*i,nper);
    end
    for j=1:Mobl_dim
        mobility_mat(j,1,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(1,1))/(npel/Mobl_dim);
        for i=2:Mobl_dim-1
            mobility_mat(j,i,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(i-1,1) ...
            & sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)<=bound_btm(i,1))/(npel/Mobl_dim);
        end
        mobility_mat(j,Mobl_dim,k)=sum(sort_mid_nper((j-1)*npel/Mobl_dim+1:j*npel/Mobl_dim,1)>bound_btm(Mobl_dim-1,1))/(npel/Mobl_dim);
    end
    Shorrocks_S(k,1) = (Mobl_dim-trace(mobility_mat(:,:,k)))/(Mobl_dim-1);
end
end

function figHandle = create_new_figure(figTitle)
persistent figCounter;
if isempty(figCounter)
    figCounter = 0;
end

figCounter = figCounter + 1;
timestamp = datestr(now,'yyyymmddTHHMMSSFFF');
uniqueName = sprintf('%s (%d - %s)', figTitle, figCounter, timestamp);
figHandle = figure('Name', uniqueName, 'NumberTitle','off');
end

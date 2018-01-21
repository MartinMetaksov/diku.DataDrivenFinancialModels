
% a) diversification
disp('a) diversification');

% Microsoft, Boeing, Bank of America, Exxon Mobil, UnitedHealth Group,
% Pepsico, NextEra Energy
tickers = {'MSFT' 'BA' 'BAC' 'XOM' 'UNH' 'PEP' 'NEE'};
date_start = '01121997';
date_end = '31122017';

% Prevent unnessary loading of data from yahoo finance
if ~exist('stocks', 'var')
    stocks = hist_stock_data(date_start, date_end, tickers, 'frequency', 'mo');
end
disp(['Collected stock data (', date_start, ' - ', date_end, ') of ', strjoin(tickers)]);

% b) estimate
disp('b) estimate');

% Caclulate log returns with a rolling window of 10 years
ymeans = cell(1,11);
ystds = cell(1,11);
ycorrs = cell(1,11);

for wi = 1:11
    
    LogReturns = zeros(120, length(stocks));
    for si = 1:length(stocks)
        s = stocks(si);
        pPrice = s.AdjClose((wi-1)*12+1:(wi+9)*12);
        cPrice = s.AdjClose((wi-1)*12+2:(wi+9)*12+1);
        LogReturns(:, si) = log(cPrice ./ pPrice);
    end
    
    ymeans{wi} = 12 * mean(LogReturns)';
    ystds{wi} = sqrt (12 * var(LogReturns))';
    ycorrs{wi} =  corr(LogReturns)';
end
disp(['Estimated expected yearly returns for ', num2str(length(ymeans)), ' windows with ', num2str(length(stocks)), ' stocks']);

% c) efficient frontier
% d) Tobin separation
disp('c) efficient frontier');
disp('d) Tobin separation');

RF = .01;
xopt = cell(11);
xopt2 = cell(11);
muopt = zeros(11,1);
muopt2 = zeros(11,1);
sigopt = zeros(11,1);
sigopt2 = zeros(11,1);
C = cell(11,1);

figure('Name', 'Efficient Frontier');
for i=1:11
    [xopt{i}, muopt(i), sigopt(i)]  = highest_slope_portfolio( ycorrs{i}, RF, ymeans{i}, ystds{i} );
    subplot(4,3,i);
    hold on
    title(i)
    xlim([0, 0.5])
    ylim([-0.2, 0.4])
    xlabel('\sigma_p')
    ylabel('$\bar{R}_p$','Interpreter','Latex')
    plot(sigopt(i), muopt(i) , 'x')
    RF_p1 = [0 sigopt(i) 2* sigopt(i)];
    opt1_p = [.01  muopt(i) (2 * muopt(i) - RF) ];
    line(RF_p1, opt1_p  );
    C{i}=diag(ystds{i})*ycorrs{i}*diag(ystds{i});
    [xopt2{i}, muopt2(i), sigopt2(i)]  = highest_slope_portfolio( ycorrs{i}, 0.02, ymeans{i}, ystds{i} );
    hold off
end

large_n = 100;
k = 20;
mu_p = zeros(4, 4*k* large_n + 1);
std_p = zeros(4, 4*k* large_n + 1);

for j=1:11
    for i = -2*k* large_n:1:2 * 2*k*large_n
        curr_port = i / large_n * xopt2{j} + (1 - i / large_n) * xopt{j};
        mu_p (j, i + 2*k * large_n + 1) = curr_port' * ymeans{j};
        std_p(j, i + 2*k * large_n + 1) = sqrt(curr_port' * C{j} * curr_port);
        
    end
    subplot(4,3,j);
    hold on
    plot( std_p(j,:), mu_p(j,:));
    hold off
end
disp(['Plotted ', num2str(length(ymeans)), ' efficient frontiers with Tobin separation and risk free rate ', num2str(RF)]);

% e) asset allocation
%c represent goal returns for the asset allocation, change this and notice
%the differences in the backtest
disp('e) asset allocation');
c = 0.15; 
disp(['Computing optimal portfolios for goal returns of ', num2str(c)]);

portf=zeros(11,8);

for i=1:11
    syms x;
    r = solve(x*RF+(1-x)*muopt(i)==c, x);
    
    for j=1:8
        if (j<=7)
            portf(i,j)=(1-r)*xopt{i}(j);
        else
            portf(i,j)= r;
        end
    end
end
    
to=zeros(1,10);
for i=1:10
    to(i)=sum(abs(portf(i,:)-portf(i+1,:))*100);
end

disp(['Average portfolio turnover is ', num2str(mean(to))]);

% f) backtest
disp('f) backtest');
backtest = zeros(10,8);
for yi = 11:20 %for year 11-20
    
    for si = 1:length(stocks)
        s = stocks(si);
        cPrice = s.AdjClose(yi*12);
        pPrice = s.AdjClose((yi-1)*12);
        backtest(yi-10, si) = portf(yi-10, si) * log(cPrice / pPrice);
    end
    backtest(yi-10, 8) = portf(yi-10, 8) * RF;
end

%Returns for each year 10-11,11-12... using the portfolios from the asset
%allocation

backtest_yearly = sum(backtest, 2);
backtest_yearly_avg = mean(backtest_yearly);
backtest_yearly_dev = sqrt(var(backtest_yearly));

% Show average return and the standard deviation 
disp(['Average portfolio return is ', num2str(backtest_yearly_avg)]);
disp(['Standard deviation is ', num2str(backtest_yearly_dev)]);

% g) beta:
% use a broad stock index to test, whether our portfolio is inline with the 
% CAPM prediction
% collect a broad stock index to test with
disp('g) beta');
if ~exist('SP', 'var')
    SP = hist_stock_data('01122007', '01012018','^GSPC','frequency','mo');    
end

SPLR = log(SP.AdjClose(13:12:end) ./ SP.AdjClose(1:12:end-12));

% Compute return on market
RM = mean(SPLR);
disp(['Average S&P market return is ', num2str(RM)]);

disp('Computing linear regression of the portfolio and S&P');
RmtmRft = SPLR - RF;
RitmRft = backtest_yearly - RF;
X = [ones(size(RmtmRft)), RmtmRft];
[b,bint,r,rint,stats] = regress(RitmRft, X); 

disp(['Betas are ', num2str(b')]);

% Calculate portfolio return using CAPM prediction
RP_test = RF + (RM - RF) * b(2);
disp(['Portfolio return using CAPM prediction is ', num2str(RP_test)]);

% Calculate Jensen's alpha
alpha = backtest_yearly_avg - RP_test;
disp(['Jensen''s alpha is ', num2str(alpha)]);

% i) timing
% Calculate Treynor-Mazuy measure
% (Rit - RFt ) = ai + bi (Rmt - RFt ) + ci (Rmt - RFt)^2 + eit
disp('i) timing');

RmtmRft2 = RmtmRft.^2;

X2 = [ones(length(RmtmRft),1), RmtmRft, RmtmRft2];

[b2,bint2,r2,rint2,stats2] = regress(RitmRft, X2);

disp(['Regression result is ', num2str(b2')]);

figure('Name', 'Timing')
scatter(SPLR, backtest_yearly)
title('Timing')
xlabel('R_m')
ylabel('R_i')
hold on
axis manual
syms x;
f1 = b(1) + b(2) * x;
f2 = b2(1) + b2(2) * x + b2(3) * x^2;
fplot(f1);
fplot(f2);
legend('Annual returns', sprintf('y = %0.2f %+0.2fx', b), sprintf('y = %0.2f %+0.2fx %+0.2fx^2', b2))
hold off


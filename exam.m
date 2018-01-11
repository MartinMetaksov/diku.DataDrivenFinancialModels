
% Prevent unnessary loading of data from yahoo finance
if ~exist('stocks', 'var')
    % Microsoft, Boeing, Bank of America, Exxon Mobil, UnitedHealth Group,
    % Pepsico, NextEra Energy
    tickers = {'MSFT' 'BA' 'BAC' 'XOM' 'UNH' 'PEP' 'NEE'};
    
    stocks = hist_stock_data('01011998', '01012018', tickers, 'frequency', 'mo');
end

% Caclulate log returns with a rolling window of 10 years
ymeans = cell(1,11);
ystds = cell(1,11);
ycorrs = cell(1,11);

for wi = 0:10
    
    LogReturns = zeros(119, length(stocks));
    for si = 1:length(stocks)
        s = stocks(si);
        LogReturns(:, si) = log(s.AdjClose(wi*12+2:wi*12+120) ./ s.AdjClose(wi*12+1:wi*12+119));
    end
    
    ymeans{wi+1} = 12 * mean(LogReturns)';
    ystds{wi+1} = sqrt (12 * var(LogReturns))';
    ycorrs{wi+1} =  corr(LogReturns)';
end

%calculate the efficient frontiers

RF = .01;
xopt = cell(11);
xopt2 = cell(11);
muopt = zeros(11,1);
muopt2 = zeros(11,1);
sigopt = zeros(11,1);
sigopt2 = zeros(11,1);
C = cell(11,1);
for i=1:11
    [xopt{i}, muopt(i), sigopt(i)]  = highest_slope_portfolio( ycorrs{i}, RF, ymeans{i}, ystds{i} );
    subplot(4,3,i);
    hold on;
    title(i)
    plot (sigopt(i), muopt(i) , 'x');
    RF_p1 = [0 sigopt(i) 2* sigopt(i)];
    opt1_p = [.01  muopt(i) (2 * muopt(i) - RF) ];
    line(RF_p1, opt1_p  );
    C{i}=diag(ystds{i})*ycorrs{i}*diag(ystds{i});
    [xopt2{i}, muopt2(i), sigopt2(i)]  = highest_slope_portfolio( ycorrs{i}, 0.02, ymeans{i}, ystds{i} );
    hold off;
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
    hold on;
    plot( std_p(j,:), mu_p(j,:));
    hold off;
end
%c represent goal returns for the asset allocation, change this and notice
%the differences in the backtest
c=0.3; 
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


for wi = 10:19 %for year 10-19
    
    LogReturn = zeros(11, length(stocks));
    for si = 1:length(stocks)
        s = stocks(si);
        LogReturn(:, si) = log(s.AdjClose(wi*12+2:wi*12+12) ./ s.AdjClose(wi*12+1:wi*12+11));
    end
    
end
backt=zeros(8,10);
ilog=LogReturn';
iport=portf';
for i= 1:10
    backt(1:7,i)=ilog(:,i).*iport(1:7,i); %portfolio times logreturns
    backt(8,i)=iport(8,i)*RF; %Risk free allocation
end

%Returns for each year 10-11,11-12... using the portfolios from the asset
%allocation
rtrns=sum(backt)




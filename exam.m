
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

c=0.1;
v=zeros(1,11);

for i=1:11
    syms x;
    v(i)=solve(x*RF+(1-x)*muopt(i)==c, x);
end

portf=zeros(11,8);
for i=1:11
    for j=1:8
        if (j<=7)
            portf(i,j)=v(i)*xopt{i}(j);
        else
            portf(i,j)=v(i);
        end
    end
end
to=zeros(1,11);
for i=1:10
    to(i)=abs(sum(portf(i,:)-portf(i+1,:)*100));
end


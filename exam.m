
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



% 2. Bonds. Below is the link for the bloomberg bonds where you can find some 
% bonds to use. Assume that the price is the clean price for a bond. As the 
% specific dates of the bond are not written (but rather the amount of years left), 
% it is OK if you assume the dates and that each bond would pay semi annual coupons. 
% Please use bonds of maturity of at least 2 years so the price will be set correctly.

% https://www.bloomberg.com/markets/rates-bonds/government-bonds/us
EUR_TO_USD_RATE = 1.22; % per 21.01.2018

% a) Use five different "real" bonds and calculate for these bonds the yield 
%    to maturity, duration and convexity.

f1 = 'BOND_NAME'; f2 = 'MATURITY'; f3 = 'COUPON'; f4 = 'PRICE'; f5 = 'YIELD'; 
% Treasury yields
b1 = struct(f1, 'GT10:GOV', f2, 10, f3, 2.25, f4, 96.48,f5, 2.66);
b2 = struct(f1, 'GT5:GOV', f2, 5, f3, 2.13, f4, 98.50,f5, 2.45);
b3 = struct(f1, 'GT2:GOV', f2, 2, f3, 1.88, f4, 99.64,f5, 2.06);
% Treasury Inflation Protected Securities (TIPS)
b4 = struct(f1, 'GTII10:GOV', f2, 10, f3, 0.50, f4, 99.20,f5, 0.58);
b5 = struct(f1, 'GTII5:GOV', f2, 5, f3, 0.13, f4, 98.89,f5, 0.39);
Bonds = [b1; b2; b3; b4; b5];

face = 100000 * EUR_TO_USD_RATE; % using the same variable for all bonds (for simplicity)
BLength = length(Bonds);
for i = 1:BLength
    % Finding yield to maturity
    b = Bonds(i);
    C = b.COUPON * face / 100;
    syms y;
    ytm_eqns = sym(zeros(b.MATURITY+1,1));
    for t = 1:(b.MATURITY*2)
        ytm_eqns(t) = C / (1 + y/2)^t;
        if t == (b.MATURITY*2)
            ytm_eqns(t+1) = (face) / (1+y/2)^t;
        end
    end
    ytm_eqn = sum(ytm_eqns) == b.PRICE;
    ytm_real = real(double(solve(ytm_eqn, y)));
    ytm_i = find(ytm_real>0);
    ytm = ytm_real(ytm_i(1));
    [Bonds(i).YTM] = ytm;
    
    % Duration and convexity
    D = 0;
    DC = 0;
    for t = 1:(b.MATURITY*2)
        D = D + (t * C)/(1 + ytm / 2)^t;
        DC = DC + (t * (t + 1) * C)/(1 + ytm / 2)^t;
        if (t == (b.MATURITY*2))
            D = D + (t * face) / (1 + ytm / 2)^t ;
            DC = DC + (t * (t + 1) * C)/(1 + ytm / 2)^t;
        end
    end
    % Slide 4 from ch21.pdf     
    D_modified = D / face;
    D_macaulay = D_modified * (1 + ytm);
    Conv = DC / face;
    [Bonds(i).DURATION] = D;
    [Bonds(i).DURATION_MODIFIED] = D_modified;
    [Bonds(i).DURATION_MACULAY] = D_macaulay;
    [Bonds(i).D_CONVEXITY] = DC;
    [Bonds(i).CONVEXITY] = Conv;
    
    fprintf('BOND_NAME = %s; MATURITY = %s; COUPON = %s; PRICE = %s; YIELD = %s\n', ...
        b.BOND_NAME, num2str(b.MATURITY), num2str(b.COUPON), num2str(b.PRICE), num2str(b.YIELD))
    fprintf('y = %s; D = %s; D_modified = %s; D_macaulay = %s; DC = %s: C = %s\n', ...
        num2str(ytm), num2str(D), num2str(D_modified), num2str(D_macaulay), num2str(DC), num2str(Conv));
    fprintf('\n');
end

% b) Calculate the duration and convexity of a portfolio of these bonds, 
%    if EUR 100.000,-- is invested in each of them
PortfolioPrice   = 100000 * EUR_TO_USD_RATE * 5;
PortfolioWeights = ones(BLength,1)/BLength;

PortfolioDuration  = PortfolioWeights' * [Bonds.DURATION_MACULAY]';
PortfolioConvexity = PortfolioWeights' * [Bonds.CONVEXITY]';


% c) Estimate the potential decline in the market value of your portfolio, 
%    if the yield increases by 150 basis points.
dY = 0.015; % 150 points
PortfolioAmounts = PortfolioPrice * PortfolioWeights ./ [Bonds.PRICE]';
m
NBP = zeros(BLength, 1);
for i = 1:BLength
    b = Bonds(i);
    C = b.COUPON * face / 100;
    New_YTM = b.YTM + dY;
    P = 0; 
    for t = 1:(b.MATURITY*2)
        P = P + C / (1 + New_YTM / 2)^t;
        if t == (b.MATURITY*2)
            P = P + (face) / (1+ New_YTM / 2)^t;
        end
    end
    NBP(i) = P;
end

NewPortfolioPrice = PortfolioAmounts' * NBP;
s = PortfolioPrice - NewPortfolioPrice;
diff = abs(s);
change = diff / PortfolioPrice * 100;
change_s = 'a decrease';
if diff < 0
    change_s = 'an increase';
end
fprintf('old portfolio price - new portfolio price = $%s, which is %s of %s%%\n', num2str(diff), change_s,num2str(change));

function [ port, opt_mu, opt_sigma ] = highest_slope_portfolio( R, RF, mu, sigma )
    % This function finds the portfolio with the largest slope
    % this function can easily be much more general 
    % e.g. mu, RF, sigma can be parameters
    if  nargin < 1
        return
    elseif nargin == 1
        RF=0.02;
        mu=[.1 .2]';
        sigma=[.1 .2]';
    elseif nargin == 2
        mu=[.1 .2]';
        sigma=[.1 .2]';
    end

    % Here we use our define correlation coefficient
    C=diag(sigma)*R*diag(sigma);

    A=2*C;
    A(:,end+1)=-(mu-RF);
    A(end+1,1:end-1)=(mu-RF)';

    Rp=mu(1);
    b=zeros(length(mu),1);
    b(end+1,1)=Rp-RF;

    x=inv(A)*b;
    xopt=x(1:length(mu))./sum(x(1:length(mu)))';

    % Return value
    port = xopt;
    opt_mu  = xopt' * mu;
    opt_sigma = sqrt( xopt' * C * xopt);

end


function dt = getRebalancingFrequency(rebalancing_frequency)
 
    trading_days_per_year = 252;
    weeks_per_year = 48;
    months_per_year = 12;
    
    if strcmpi(rebalancing_frequency, 'daily')
        dt = 1 / trading_days_per_year;
    elseif strcmpi(rebalancing_frequency, 'weekly')
        dt = 1 / weeks_per_year;
    elseif strcmpi(rebalancing_frequency, 'monthly')
        dt = 1 / months_per_year;
    else
        error('Invalid rebalancing frequency. Choose "daily", "weekly", or "monthly".');
    end
end
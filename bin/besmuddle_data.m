%% Function to add artificial baseline to data set
function dat = besmuddle_data(option, data1) 
    switch(option)
        case 1
            f  = @(y) (-0.0005*(y-5)^2 + 0.01); % 1 Percent max 
        case 2
            f = @(y) (0.0002*(y)^2 + 0.002*y+0.01); % 5 Percent max
        case 3
            f = @(y) (0.002*(y-8)^2+0.002*y+0.03); % 10 Percent max
    end
    baseline = zeros(size(data1.X{1}));
    for j=1:length(baseline)
       x = data1.X{1}(j);
       baseline(j) = f(x);
    end
    for i=1:length(data1.D)
        data1.D{i} = data1.D{i}+baseline;
        data1.ma(i) = max(data1.D{i});
    end
    data1.scal = data1.scal*max(data1.ma);
    dat = data1;
end
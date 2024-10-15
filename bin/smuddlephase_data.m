% Adds incorrect Phasing retroactively to the data
function dat = smuddlephase_data(option, data1)
    switch(option)
        case 1
            phi = pi/18; % Small amounts of shifting
        case 2
            phi = pi/9; % Medium amounts of shifting
        case 3
            phi = pi/6; % Copious amounts of shifting
    end

    for j=1:length(data1.D)
        dold        = hilbert(data1.D{j});  % Constructs imaginary part
        dfour       = ifft(dold);           % Inverse Fourier transform 
        dfour       = dfour*exp(1i*phi);    % Shifting
        dnew        = fft(dfour);           % Fourier transform 
        data1.D{j}  = real(dnew);           % Only need real part
        data1.ma(j) = max(data1.D{j});      % Adjust maxima
    end
    data1.scal = data1.scal*max(data1.ma);
    dat = data1;
end
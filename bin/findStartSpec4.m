% Overhead function for 'findMinima3', and finding peaks. Does plots as
% well.
function center = findStartSpec4(Xinit, Dinit, bw, epsi, delta, d, N, info)
    
    if size(Dinit, 1) == 1 % 'mysgfilt' only accepts column vectors
        Dinit = Dinit';
    end
    [d_der0, d_der1, d_der2, ~] = mysgfilt(d,N,Dinit); % Smoothing
    [center, d2center] = findMinima3(Dinit, d_der1, d_der2, bw, epsi, delta); % Again, third time's the charm. (Actual calculations)
    
    % Plots
    if info >= 2
       figure(1);
       clf;
       if info>=3
           sgtitle('Unsmoothed version')
       end
       subplot(1,2,1)
       title('Unsmoothed data')
       plot(Xinit, Dinit)
       hold on
       plot(Xinit(center),Dinit(center), 'r*')
       hold on
       line([Xinit(1) Xinit(end)], [epsi, epsi], 'Color', 'black');
       set(gca, 'XDir', 'reverse');
       axis tight
       
       subplot(1,2,2)
       title('Second derivatives')
       plot(Xinit, d_der2)
       hold on
       
       plot(Xinit(d2center), d_der2(d2center), 'g*')
       hold on
       line([Xinit(1), Xinit(end)], [delta,delta], 'Color', 'black');
       set(gca, 'XDir', 'reverse');
       axis tight
       
       % More plots
       if info >= 3
           figure(2);
           clf;
           sgtitle('Smoothed version')
           subplot(1,2,1)
           title('Smoothed data')
           plot(Xinit, d_der0)
           hold on
           plot(Xinit(center),d_der0(center), 'r*')
           hold on
           line([Xinit(1) Xinit(end)], [epsi, epsi], 'Color', 'black');
           set(gca, 'XDir', 'reverse');
           axis tight
           
           subplot(1,2,2)
           title('Second derivatives')
           plot(Xinit, d_der2)
           hold on
       
           plot(Xinit(d2center), d_der2(d2center), 'g*')
           hold on
           line([Xinit(1), Xinit(end)], [delta,delta], 'Color', 'black');
           set(gca, 'XDir', 'reverse');
           axis tight
       end
    end
end
function [] = plotXY(trial, angles)

    figure

    for j = 1:length(angles)
        subplot(3,3,j)
        for i = 1:100
            t = 1:length(trial(i,j).handPos(1,:));
            plot3(t,trial(i,j).handPos(1,:),trial(i,j).handPos(2,:));
            hold on
        end
        grid on
        title(['Movement Direction =' num2str(angles(j)*180/pi)])
        xlabel('time')
        ylabel('x')
        zlabel('y')
        az = 90;
        el = 0;
        view(az, el);
        ylim([-120 120])
        zlim([-120 120])
    end
end
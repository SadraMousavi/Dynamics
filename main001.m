
designParam.R = 0.25;
designParam.g = 9.81;
designParam.P = 10;
designParam.l = 2;
designParam.m = 2;

GenCoord_0.x     = 2*designParam.R;
GenCoord_0.theta = 2*atan(designParam.R/GenCoord_0.x);
GenCoord_0.phi   = -pi/2;

plotter(designParam, GenCoord_0)

[T,X] = ode45(@(t,x)Dynamics(t,x,designParam),0:0.01:0.7,[GenCoord_0.theta 0]);
GenCoord = getGenCoordFromSolution(X, designParam, GenCoord_0);
animator001(designParam, GenCoord, T)



function [dx] = Dynamics(t,x,designParam)
    dx=zeros(2,1);
    R=designParam.R;
    g=designParam.g;
    P=designParam.P;
    l=designParam.l;
    m=designParam.m;
    %
    dx(1) = x(2);
    dx(2) = - P*l - m*g*l*cos(x(1))/2 + m*R^2*x(2)^2*sin(x(2)) * ( 2+cos(x(1)) )/( (1-cos(x(1)))^3 );
    dx(2) = dx(2) / (m*l^2/3 + R^2*( (2+cos(x(1)))^2 )/( (1-cos(x(1)))^2 ) ) / m;
end

function [] = plotter(designParam, GenCoord)
    bar_x = designParam.l*cos(GenCoord.theta);
    bar_y = designParam.l*sin(GenCoord.theta);
    plot([0, bar_x], [0, bar_y])
    hold on
    %
    circle_x = GenCoord.x;
    circle_y = designParam.R;
    line1_x = circle_x + designParam.R*cos(GenCoord.phi);
    line1_y = circle_y + designParam.R*sin(GenCoord.phi);
    plot([circle_x, line1_x], [circle_y, line1_y])
    %
    theta_ = 0:pi/50:2*pi;
    x_ = circle_x + designParam.R*cos(theta_);
    y_ = circle_y + designParam.R*sin(theta_);
%     viscircles([circle_x,circle_y],designParam.R);
    plot(x_, y_)
    xlim([-2 4]);
    ylim([-2 4]);
%     axis equal
end

function [GenCoord] = getGenCoordFromSolution(x, designParam, GenCoord_0)
    GenCoord(length(x)) = struct('theta',[],'phi',[],'x',[]);
    for i=1:length(x)
        GenCoord(i).theta = x(i) ;
        GenCoord(i).phi = x(i) + cot(x(i)/2) - 1.3565 - pi;
        GenCoord(i).x = designParam.R * cot(x(i)/2);
    end
end



function [F] = animator001(designParam, GenCoord, t)
%
figure(100)
plot(0,0,'bs');
ax = gca;
ax.NextPlot = 'replaceChildren';
%
loops = length(GenCoord);
F(loops) = struct('cdata',[],'colormap',[]);
%
% hold on
for i=1:loops
    plotter(designParam, GenCoord(i))
    % text(-0.2,0, [num2str(t(i)) ' s'],'Color','r', 'FontSize', 30)
    title(sprintf('t = %1.3f s', t(i)),'Color','r', 'FontSize', 30)
    %
    hold off
    xlabel('x (m)');
    ylabel('y (m)');
    set(gca, 'fontsize', 30)
    drawnow
    F(i) = getframe;
    im{i} = frame2im(F(i));
end
% movie(F)



filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:loops
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename, 'gif','LoopCount',Inf);%'DelayTime',1
    else
        imwrite(A,map,filename,'gif','WriteMode','append');
    end
end



end






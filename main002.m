
designParam.R = 0.30;
designParam.g = 9.81;
designParam.r = 0.10;
designParam.h = 3;
designParam.m = 1;
designParam.M = 0.2;
designParam.alpha = 30*pi/180;
designParam.beta  = 40*pi/180;

GenCoord_0.x     = 6*designParam.R;
GenCoord_0.y     = 6*designParam.R;
GenCoord_0.theta = 0;
GenCoord_0.phi   = 0;

plotter(designParam, GenCoord_0)

% [T,X] = ode45(@(t,x)Dynamics(t,x,designParam),0:0.01:0.7,[GenCoord_0.theta 0]);
T = 0:0.01:0.8;
a = designParam.g * (designParam.m*sin(designParam.beta) - 2*sin(designParam.beta)*designParam.R*designParam.M/(designParam.R-designParam.r)) / (designParam.m+4*designParam.R^2*designParam.M/(designParam.R-designParam.r)^2) / 1.5;
X = 0.5*a*T.^2;
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
    edge_1 = - designParam.h/tan(designParam.alpha);
    edge_2 = + designParam.h/tan(designParam.beta);
    plot([edge_1, 0, edge_2], [0, designParam.h, 0], 'k')
    hold on
    plot([edge_1, edge_2], [0, 0], 'k')
    %
    circle1_x = - GenCoord.y*cos(designParam.alpha) - designParam.R*sin(designParam.alpha);
    circle1_y = - GenCoord.y*sin(designParam.alpha) + designParam.R*cos(designParam.alpha) + designParam.h;
    line1_x = circle1_x + designParam.R*cos(GenCoord.phi);
    line1_y = circle1_y + designParam.R*sin(GenCoord.phi);
    plot([circle1_x, line1_x], [circle1_y, line1_y], 'r')
    %
    circle2_x = + GenCoord.x*cos(designParam.beta) + designParam.r*sin(designParam.beta);
    circle2_y = - GenCoord.x*sin(designParam.beta) + designParam.r*cos(designParam.beta) + designParam.h;
    line2_x = circle2_x + designParam.r*cos(GenCoord.theta);
    line2_y = circle2_y + designParam.r*sin(GenCoord.theta);
    plot([circle2_x, line2_x], [circle2_y, line2_y], 'r')
    %
    theta_ = 0:pi/50:2*pi;
    x1_ = circle1_x + designParam.R*cos(theta_);
    y1_ = circle1_y + designParam.R*sin(theta_);
    plot(x1_, y1_, 'b')
    x2_ = circle1_x + designParam.r*cos(theta_);
    y2_ = circle1_y + designParam.r*sin(theta_);
    plot(x2_, y2_, 'b')
    x3_ = circle2_x + designParam.r*cos(theta_);
    y3_ = circle2_y + designParam.r*sin(theta_);
    plot(x3_, y3_, 'b')
    %
    point_y = 2*designParam.r/sin(pi/2-designParam.beta) + designParam.h;
    circle1_x1 = - GenCoord.y*cos(designParam.alpha) - (designParam.R-designParam.r)*sin(designParam.alpha);
    circle1_y1 = - GenCoord.y*sin(designParam.alpha) + (designParam.R-designParam.r)*cos(designParam.alpha) + designParam.h;
    circle2_x1 = + GenCoord.x*cos(designParam.beta) + 2*designParam.r*sin(designParam.beta);
    circle2_y1 = - GenCoord.x*sin(designParam.beta) + 2*designParam.r*cos(designParam.beta) + designParam.h;
    plot([circle1_x1, 0, circle2_x1], [circle1_y1, point_y, circle2_y1])
    plot([0, 0], [designParam.h, point_y], 'k')
    plot([0, 0], [designParam.h, point_y], 'k')
    plot(0, point_y, 'ok')
    %
    xlim([-6 4]);
    ylim([-3 7]);
%     axis equal
end

function [GenCoord] = getGenCoordFromSolution(x, designParam, GenCoord_0)
    GenCoord(length(x)) = struct('theta',[],'phi',[],'x',[],'y',[]);
    for i=1:length(x)
        GenCoord(i).theta = -x(i) / designParam.r + GenCoord_0.theta;
        GenCoord(i).phi = 2*designParam.R/(designParam.R-designParam.r) * GenCoord(i).theta + GenCoord_0.phi;
        GenCoord(i).x = x(i) + GenCoord_0.x;
        GenCoord(i).y = designParam.R * GenCoord(i).phi + GenCoord_0.y;
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






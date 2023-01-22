Preload = 0:0.0001:15;
N = size(Preload);
N = N(2);
ControlLine = (10.3+(-10.3./(1+(Preload/7).^2.3)));
plot(Preload,ControlLine,'b','linewidth',2)
hold on
Y1 = -1.96.*(Preload - 6.774) + 4.955;
Y2 = -1.96.*(Preload - 10) + 4.955;
Y3 = -1.96.*(Preload - 3.5) + 4.955;
plot(Preload,Y1,'r','linewidth',2)
plot(Preload,Y2,'r','linewidth',2)
plot(Preload,Y3,'r','linewidth',2)
plot(6.774,4.955,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
plot(10,4.955,'o-','MarkerFaceColor','red','MarkerEdgeColor','black')
plot(9.124,6.672,'o-','MarkerFaceColor','green','MarkerEdgeColor','black')
plot(3.5,4.955,'o-','MarkerFaceColor','red','MarkerEdgeColor','black')
plot(4.586,2.826,'o-','MarkerFaceColor','green','MarkerEdgeColor','black')
axis([0 15 0 10])
xlabel('Ventricular Preload (mmHg)')
ylabel('Pump Flow (L/min)')

x1 = [6.95 9.78];
y1 = [4.955 4.955];
x2 = [9.95 9.17];
y2 = [5.1 6.55];
x3 = [6.6 3.7];
y3 = [4.955 4.955];
x4 = [3.55 4.52];
y4 = [4.85 2.95];

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
drawArrow(x1,y1,'linewidth',3,'color','k')
drawArrow(x2,y2,'linewidth',3,'color','k')
drawArrow(x3,y3,'linewidth',3,'color','k')
drawArrow(x4,y4,'linewidth',3,'color','k')

% CRC = zeros(1,N);
% 
% for i=1:N
%     if Preload(i) > 0
%         CRC(i) = (10.5+(-10.5./(1+(Preload(i)/4).^2.7)));
%     end
% end
% 
% % Y = (10.5+(10.5./(1+(Preload/4).^2.7)));
% Y = sigmoide(Preload,0);
% x1 = 5;
% y1 = 3;
% slope = -0.65;
% y = slope.*(Preload - x1) + y1;
% % epsilon = 1e-6;
% % idx = find(y - CRC < epsilon, 1); %// Index of coordinate in array
% % x_target = Preload(idx);
% % y_target = y(idx);
% plot(Preload,CRC,Preload,Y)
% axis([-4 14 0 12])
% 
% function f = sigmoide(x,A)
%     for i = 1:length(x)
%         if x(i) < A
%             f(i) = 6.25;
%         elseif x(i) >= A
%             f(i) = -0.65.*(x(i) - 5) + 3;
%         end
%     end
% end

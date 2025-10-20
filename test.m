k1 = 10; k2=1;
x = 0:0.001:1;
p = 2;
q = 0.5;
r = -0.5;
y1 = x*k1+(1-x)*k2;
y2 = (k1^q*x+k2^q*(1-x)).^(1/q);
y3 = (k1^r*x+k2^r*(1-x)).^(1/r);
y4 = 1./((1/k1-1/k2)*x+1/k2);
y5 = x.^3*k1+(1-x.^3)*k2;
plot(x,y1,'r', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on
plot(x,y2,'g', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on
plot(x,y3,'b', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on
plot(x,y4,'c', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on
plot(x,y5,'k', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('p = 1','p=0.5','p=-0.5','p=-1','SIMP:3');
title('系数族曲线对比', ...
      'FontSize', 16, 'FontWeight', 'bold');
% 设置网格
grid on;
grid minor;

% 优化坐标轴
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;

text(0.02, 0.98, ...
     {sprintf('k1 = 10, k2 = 1')}, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', ...
     'FontSize', 10);
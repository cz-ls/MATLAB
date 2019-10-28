% ������λ�ã����ڶԳ��ԣ������ڵڶ����޵�1/4Բ���ϣ�
m = 100001;

n1 = 1;
n2 = 1.3;

theta = linspace(pi/2,pi,m); 

% ���ؿ��޷� theta = rand(1,m)*pi/2+pi/2;

F = zeros(180,1);

for i = 2:m-1
    
%     % ���ؿ��޷�
%     if theta(i) == pi/2 || theta(i) == pi
%         continue;
%     end
    
    theta_i = pi-theta(i);
    theta_t = asin(sin(theta_i)*n1/n2);

    rs = ((n1*cos(theta_i)-n2*cos(theta_t))/(n1*cos(theta_i)+n2*cos(theta_t)))^2;
    rp = ((n1*cos(theta_t)-n2*cos(theta_i))/(n1*cos(theta_t)+n2*cos(theta_i)))^2;
    R = (rs+rp)/2;
    rest = 1;

    % ��¼��ʼ�����ĽǶȡ�����
    k = ceil((pi-theta_i*2)/pi*180);
    F(k) = F(k)+1*R;

    % ��·��ͼ�α�Ե�Ӵ���ĸ���n�����������ڳ�ʼֵ��0.1%����׷�٣�
    n = 0;
    while rest>0.001
        if n == 0
            rest = rest*(1-R);
        else
            rest = rest*R;
        end
        n = n+1;
    %     disp(rest);
    end

    spot = zeros(n,2);  % ��¼�Ӵ����λ��
    spot(1,:) = [cos(theta(i)) sin(theta(i))];

    outline = zeros(n,2); % ��¼������ߵ�λ��
    outline(1,:) = [spot(1,1)+cos(pi-theta_i*2) spot(1,2)+sin(pi-theta_i*2)];

    for j = 2:n
        normal_theta = atan(abs((spot(j-1,2))/(spot(j-1,1))));
        if spot(j-1,1)>0
            if spot(j-1,2)<0
                normal_theta = 2*pi-normal_theta;
            end
        else
            if spot(j-1,2)>0
                normal_theta = pi-normal_theta;
            else
                normal_theta = pi+normal_theta;
            end
        end

        temp_theta = normal_theta-(pi-2*theta_t)/2;
        temp_l = sin(theta_t);
        temp = [temp_l*cos(temp_theta) temp_l*sin(temp_theta)];
        spot(j,:) = [2*temp(1)-spot(j-1,1) 2*temp(2)-spot(j-1,2)];

        % ��¼�����������ĽǶȡ�����
        if j>2
            location_scatter = normal_theta-theta_i;
            outline(j-1,:) = [spot(j-1,1)+cos(location_scatter) spot(j-1,2)+sin(location_scatter)];
            theta_scatter = acos(cos(location_scatter))/pi*180;
%             disp(theta_scatter);
            k = ceil(theta_scatter);
            F(k) = F(k)+(1-R)*R^(j-2);
        end
    end

    % ��¼���һ�������ĽǶȡ�����
    j = n+1;
    normal_theta = atan(abs((spot(j-1,2))/(spot(j-1,1))));
    if spot(j-1,1)>0
        if spot(j-1,2)<0
            normal_theta = 2*pi-normal_theta;
        end
    else
        if spot(j-1,2)>0
            normal_theta = pi-normal_theta;
        else
            normal_theta = pi+normal_theta;
        end
    end
    location_scatter = normal_theta-theta_i;
    outline(j-1,:) = [spot(j-1,1)+cos(location_scatter) spot(j-1,2)+sin(location_scatter)];
    theta_scatter = acos(cos(location_scatter))/pi*180;
%     disp(theta_scatter);
    k = ceil(theta_scatter);
    F(k) = F(k)+(1-R)*R^(j-2);

end

fic = F/sum(F);
figure;
plot((1:180),fic);

set(gca, 'Fontname', 'Times New Roman','FontSize',14);

xlabel('Scatter Angle [��]','Fontname', 'Times New Roman','Fontsize',15);
ylabel('Phase Function','Fontname', 'Times New Roman','Fontsize',15);

h = title(['n2 = ' num2str(n2)]);
set(h,'Fontname', 'Times New Roman','Fontsize',15);

set(gcf,'position',[100,100,700,500]);
set(gcf,'paperpositionmode','auto');
% print -dpng -r600 hw01;
print -dpng -r600 test001;

% % ���������

% % ͼ�εĶ��塢ͼ��
% r = 1;
% omega = linspace(0,2*pi);
% x = r*cos(omega);
% y = r*sin(omega);
% 
% figure;
% plot(x,y,'LineWidth',2,'Color',[0,0.4471,0.7412]);
% axis equal;
% hold on;
% grid on;

% plot([spot(1,1)-1 spot(1,1)],[spot(1,2) spot(1,2)],'-',...
%     'Color',[0.9255,0.6627,0.5490],...
%     'LineWidth',2);
% hold on;
% 
% % ͼ���ڵ�һ����
% plot([spot(1,1) spot(2,1)],[spot(1,2) spot(2,2)],'-',...
%     'Color',[0.9255,0.6627,0.5490],...
%     'LineWidth',2*(1-R));
% hold on;
% % ��һ�����䷴��⣨ɢ�䣩
% plot([spot(1,1) outline(1,1)],[spot(1,2) outline(1,2)],'-',...
%     'Color',[0.9255,0.6627,0.5490],...
%     'LineWidth',2*R);
% hold on;
% 
% for z = 2:n-1   % ͼ���ڷ����
%     plot([spot(z,1) spot(z+1,1)],[spot(z,2) spot(z+1,2)],'-',...
%         'Color',[0.9255,0.6627,0.5490],...
%     'LineWidth',2*(1-R)*R^(z-1));
%     hold on;
% end
% for z = 2:n     % ��������⣨ɢ�䣩
%     plot([spot(z,1) outline(z,1)],[spot(z,2) outline(z,2)],'-',...
%         'Color',[0.9255,0.6627,0.5490],...
%         'LineWidth',2*(1-R)*R^(z-2)*(1-R));
%     hold on;
% end
% for z = 1:n     % ����
%     plot([0 spot(z,1)*2],[0 spot(z,2)*2],'--',...
%         'Color',[0,0.4471,0.7412],...
%         'LineWidth',1);
%     hold on;
% end
% 
% plot(0,0,'.','MarkerSize',15);
% hold on;
% plot(spot(:,1),spot(:,2),'.','MarkerSize',15,'Color',[0.8510,0.3255,0.0980]);
% hold on;
% 
% set(gcf,'position',[100,100,900,900]);
% set(gcf,'paperpositionmode','auto');
% print -dpng -r300 test01;
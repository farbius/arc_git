function ShowMeTheTrace(Results, traces, i)
% Отображает на графике обнаруженную трассу с заданным номером

global Nwindow

figure(999); whitebg('k');

% Выбор цвета очередной трассы
color = jet; %rand(1,3);
color = color(mod(8*i,57)+7,:);

% Показываем формуляр "сигнала"
x_point = traces(1,i).Res_i;
y_point = Results(traces(1,i).Res_i).Begs(traces(1,i).Num_in_Res);

if (y_point > Nwindow/4),
    alignm = 'Top';
    y_point = y_point;%-10;
else
    alignm = 'Bottom';
    y_point = y_point;%+10;
end
text('Position',[x_point y_point],...
    'String', sprintf('%g', i),...
    'BackgroundColor','Yellow',...
    'EdgeColor',color,...
    'Color','black',...
    'Interpreter','tex','VerticalAlignment',alignm);

% Прорисовываем трассу
l = 1;
while (traces(l,i).Res_i > 0),
    k        = traces(l,i).Res_i;
    tmp      = traces(l,i).Num_in_Res;
    center_k = round(traces(l,i).Center);
    %wdth   = Results(k).Widths(tmp);
    %plot([k k], [Results(k).Begs(tmp) Results(k).Fins(tmp)],'y','LineWidth',2);
    if (traces(l,i).Clutter == 1),
        h = plot(k, center_k,'yo','LineWidth',2);
    else
        h = plot(k, center_k,'y*','LineWidth',2);
    end
    set(h,'Color',color)
    l = l+1;
end



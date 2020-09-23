%Организация          : ?
%Дата созд.           : 08.09.2015
%Дата правки          : 22.04.2016
%Версия               : 1.2	    
%Наименование         : Функция, имитирующая часть алгоритма,
%                       предназначенную для переноса в процессор
%Комментарии          : Часть I должна выполнятся каждый раз по получении 
%                       нового пакета данных от FPGA и должна выполнятся 
%                       не более 41мкс. 
function [Results] = Soft_part_I(Sf,Begs,Fins,Width)

isMlab = 1;
Nlength = length(Sf);
% ===== Параметры FPGA-алгоритма =========================================
% Подготовка формуляров
if (any(Begs) == 0)  % Превышения обнаружены?
    % Если нет, то заполняем пустой формуляр для очередного окна
    Results.Begs    = [];
    Results.Fins    = [];
    Results.Types   = [];
    Results.Centers = [];
    Results.Widths  = [];
    Results.Num     = 0;
    Results.n_iter  = 0;
    return % Выходим из функции
end

% ==== Предварительная обработка формуляров (объединение превышений) =
n_iter = 0; % Процедура итерационная!
while 1
    Lens = Fins-Begs; % Длины отрезков превышений
    if (length(Lens) == 1)
        break; % Если обнаружен единственный интервал, то выходим из while
    end
    Bets   = zeros(length(Lens)-1,1); % Вектор расстояний между отрезками превышений
    IsJoin = zeros(length(Lens)-1,1); % Флаг устранения промежутков
    for i = 1:length(Lens)-1
        %         Bets(i) = Begs(i+1)-Fins(i); % Длины промежутков
        %         % Критерий объединения:
        %         % если расстояние между соседними отрезками меньше ширины
        %         % одного из отрезков, то пара отрезков объединяется в один
        %         if ((Begs(i+1) > 0) && (Fins(i) > 0))
        %             if ((Lens(i) >= Bets(i))||(Bets(i) <= Lens(i+1))),
        %                 IsJoin(i) = 1; % Помечаем промежуток на удаление
        %             end
        %         end
        
        Spare = 5.0; % Что бы не "заливал" уж очень большие промежутки
        Bets(i) = Spare*(Begs(i+1)-Fins(i)); % Длины промежутков
        % если расстояние между соседними отрезками меньше ширины
        % обоих отрезков с запасом, то пара отрезков объединяется в один
        if ((Begs(i+1) > 0) && (Fins(i) > 0))
            if ((Lens(i) >= Bets(i)) || (Bets(i) <= Lens(i+1)))
                IsJoin(i) = 1; % Помечаем промежуток на удаление
            end
        end
    end
    
    if (sum(IsJoin) == 0)
        break; % Если больше нечего объединять, то выходим из while
    end
    
    % Производим объединение отрезков
    I = find(IsJoin == 1);
    Begs(I+1) = []; % Удаляем фронты с номерами I+1
    Fins(I)   = []; % Удаляем спады с номерами I
    
    n_iter = n_iter+1;
    if (n_iter > 50)
        % Слишком много превышений, похоже - ШПС, но пока просто выходим
        warning('The number of iterations exceeds 20!');
        break;
    end
end % Уходим на следующую итерацию

% Оценка средней частоты и ширины очередного промежутка
Centers = zeros(size(Begs));
Widths  = zeros(size(Begs));
for i = 1:length(Begs)
    x = Begs(i):Fins(i);
    Centers(i) = sum(x'.*Sf(x).^2)/sum(Sf(x).^2);
    %Widths(i)  = sum((x'-Centers(i)).^2.*Sf(x).^2)/sum(Sf(x).^2);
    Widths(i)  = length(x);
    if (Widths(i) > Nlength)
        Widths(i) = Nlength; % Идиотская формула для СКО ширины спектра!
    end
end

% Классификация по широкополосности
% I = ((Fins-Begs) > Width);
I = (Widths > Width);

% Заполняем формуляр для текущего окна
Results.Begs    = Begs;
Results.Fins    = Fins;
Results.Types   = I;
Results.Centers = Centers;
Results.Widths  = Widths;
Results.Num     = length(Begs);
Results.n_iter  = n_iter;

% if (length(Results(k+isMlab).Begs) > 1) || (length(Results(k+isMlab).Fins) > 1),
%     error('Length of Begs or Fins more than 1!');
% end

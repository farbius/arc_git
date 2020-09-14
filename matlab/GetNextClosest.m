function [i2,n2,isClutter] = GetNextClosest(Results1,k,Results2)
% Функция осуществляющая поиск ближайщего соседа для текущей отметки 
% (Results1), в следующем окне(Results2(1)), или через одно (Results2(2))
% окне по заданному критерию и допуску deltas
% Если обнаруживается, что длина найденной продолжающей отметки в
% clutterTrh раз или более превосходит длину текущей отметки, то взводится
% признак обнаружения импульсной помехи (isClutter = 1).

%
% Возвращает:
%   i2  - смещение найденной отметки (0 - не найдено продолжение)
%   n2  - номер превышения в спектре (от 1 до Results2(i2).Num)

global deltas clutterTrh

i2 = 0; n2 = 0; % Изначально ничего не найдено

isClutter = 0; % Помех тоже пока не видно

found = []; % Сюда складываем кандидатов

Begn = Results1.Begs(k); % Начало текущей отметки
Finh = Results1.Fins(k); % Конец текущей отметки
Wdth = Results1.Widths(k); % Ширина текущей отметки

% Сначала исследуем на один шаг вперёд
for k1 = 1:Results2(1).Num,
    if (Results2(1).isInTrace(k1) == 0),
        Begn1 = Results2(1).Begs(k1); % Начало текущей отметки
        Finh1 = Results2(1).Fins(k1); % Конец текущей отметки
        % Отметки перекрываются с учётом допуска (deltas(1))?
        if (Begn*(1-deltas(1)) > Finh1) || (Finh*(1+deltas(1)) < Begn1),
            % Текущая отметка с номером k1 не проходит - переходим к следующей
        else
            % Есть чем продолжить - заносим отметку в found
            found = [found; k1];
        end
    end
end

if ~isempty(found),
    i2 = 1; 
    if length(found) == 1,
        n2 = found; % Продолжить можно только одной отметкой - всё просто.
    else
        % Можно продолжить несколькими отметками, тогда всё сложнее:
        % выбираем ближайшую по длине к медианной длине трассы
        %[~,n2] = min(abs(Results2(1).Widths - med_len);
        % выбираем ближайшую по длине к длине предыдущей
        %[~,n2] = min(abs(Results2(1).Widths - Results1.Widths(k))); -- не
        %прошло проверку!
        [~,n2] = min(abs(Results2(1).Centers - Results1.Centers(k)));
    end
    if (Results2(1).Widths(n2)/Wdth >= clutterTrh),
        isClutter = 1; % Есть подозрение, что продолжаем помехой
    end
    return % Поиск закончен
end

% Если не найдено ничего среди соседей - смотрим на два шага вперёд
for k1 = 1:Results2(2).Num,
    if (Results2(2).isInTrace(k1) == 0),
        Begn1 = Results2(2).Begs(k1); % Начало текущей отметки
        Finh1 = Results2(2).Fins(k1); % Конец текущей отметки
        % Отметки перекрываются с учётом допуска (deltas(2))?
        if (Begn*(1-deltas(2)) > Finh1) || (Finh*(1+deltas(2)) < Begn1),
            % Текущая отметка с номером k1 не проходит - переходим к следующей
        else
            % Есть чем продолжить - заносим отметку в found
            found = [found; k1];
        end
    end
end

if ~isempty(found),
    i2 = 2;
    if length(found) == 1,
        n2 = found;
    else
        % Можно продолжить несколькими отметками,тогда -
        % выбираем ближайшую по расположению (Center)
        % По идее, кстати, дОлжно помогать проскакивать "помехи"
        [~,n2] = min(abs(Results2(2).Centers - Results1.Centers(k)));
    end
    if (Results2(2).Widths(n2)/Wdth >= clutterTrh),
        isClutter = 1; % Есть подозрение, что продолжаем помехой
    end
end

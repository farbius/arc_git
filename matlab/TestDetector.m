% Main, v.4.0
clc
clear
close all

%%
% FLAGS
ONE_SPECTRUM   =  0;
FULL_BANDWIDTH =  0;
ADD_ADC_RAW    =  0;
INPUT_GEN      =  1;

LOAD_MAT       =  1;

isMlab = 1;
global Fs deltas Nwindow clutterTrh Nstep
% ========================================================================
% Начальные постоянные 
% ========================================================================
Fs          = 800e6; %| Частота оцифровки
Nwindow     = 8192;  %| Величина окна БПФ
Width       = 40;    %| Порог принятия решения о "широкополосности" (в точках БПФ)
TrhManual   = 1;     %| Ручная установка порога (1 - да, 0 - автомат)
MaxMeanPlts = 7.0;   %| Максимальное среднее кол-во отдельных превышений
delta1      = 0.25;  %| Максимальное разрешённое смещение отметок на 1 шаг
delta2      = 0.5;   %| Максимальное разрешённое смещение отметок на 2 шага
Ncptrd      = 30;    %| Минимальное кол-во отметок для завязки трассы
MaxTracks   = 30;    %| Максимальное кол-во трасс (во избежания перегрузки)
MaxTraces   = 5;     %| Максимальное кол-во трасс при распозновании
clutterTrh  = 3.0;   %| См. процедуру GetNextClosest
Nlin = 3; %| Кол-во отметок, по которым корректируется помеховая отметка
% ========================================================================
% Примечание: вообще говоря, delta1 и delta2 должны вычисляться
% теоретически по заданной скорости модуляции и входному ОСШ. Кстати, это
% позволяет оптимально принимать решение о наличии ШП-сигнала при
% обнаружении "россыпи" узкополосных отметок (это сейчас делает процедура
% Soft_part_I, но по эвристическому алгоритму)
fileID = fopen('fft_6_dma.bin');
A = fread(fileID, 'uint32');
fclose(fileID);
C = zeros(length(A), 1);

f = (1:4096).*400e6 / 4096;

if INPUT_GEN == 0
CC = dec2hex(A);
if LOAD_MAT == 0

if FULL_BANDWIDTH == 1
FFT_ARRAY = zeros(4096, 1);
else
FFT_ARRAY = NaN(4096, 1);    
end

n = 0;
for i = 1 : length(A)
   C(i) =  bitand(A(i), 16383);
   
   if( C(i) == 0)
   n = n + 1; 
   else
       
   if(ONE_SPECTRUM == 1)
       if n == 1
           bin_num = bitand(A(i), 16383, 'uint32');
          if(bin_num ~= 16383)
                FFT_ARRAY(bin_num) = bitshift(A(i), -14, 'uint32'); 
          end
       end % n
   else
        bin_num = bitand(A(i), 16383, 'uint32');
          if(bin_num ~= 16383)
                FFT_ARRAY(bin_num) = bitshift(A(i), -14, 'uint32'); 
                
          end
   end % ONE_SPECTRUM
   end
end % for

% Если не технологический режим - прокатываем весь фрагмент
Nstep = n;   
n = 0;
Y = zeros(4096, Nstep);
for i = 1 : length(A)
   C(i) =  bitand(A(i), 16383);
   
   if( C(i) == 0)
   n = n + 1; 
   else
        bin_num = bitand(A(i), 16383, 'uint32');
          if(bin_num ~= 16383)
                Y(bin_num, n) = bitshift(A(i), -14, 'uint32');  
          end
   end
end % for

save('Sf.mat','Y');

else
    
load('Sf.mat');

Nstep = size(Y,2);


end % LOAD MAT


else % INPUT_GEN
%%
% signal generator

    tau     = 245.76e-6;%245.76e-6;   % period
    F2      = 74e6; % 237e6;    % start frequency
    F1      = 73.275e6;     % stop frequency 
    Ns      = ceil(tau*Fs);
    ts      = (0:Ns-1)./Fs;
    
% quadratic linear
    s       = chirp(ts,F1,tau,F2,'linear'); % convex concave ,[],'concave'
           
%    pspectrum(s,Fs,'spectrogram','FrequencyLimits',[0 400e6],'TimeResolution',1e-6);
 
    s = repmat(s, [1 10]);
   
%     snr = -20;
%     snr = 10^(-snr/20);
%       n = randn(length(s),1);
%      nn = n + std(n);
%       s = s + snr*nn;


    Nstep = ceil(length(s)/ Nwindow);
    Y = zeros(Nwindow/2, Nstep);
    for i = 0:Nstep-2
        
        sF        = abs(fft(s(i*Nwindow + 1 : (i+ 1)*Nwindow)));
        Y(:, i+1) = sF(1:Nwindow/2);
        
    end  
    
%     Nstep = 8192;
%     Y = Y(:, 17:Nstep+16);
end % INPUT_GEN



figure(999)
imagesc(Y)
title('frequency-time-power')
colorbar
hold on

% 
% figure(10)
% plot(1:4096, 20*log10(Y(:, 1:4)));


% Формуляры для каждого окна с результатами первышений,
% объединений и классификации
% Поля: начала отрезков с превышеними, концы -//-, тип, центры отрезков,
% общее кол-во отрезков, кол-во итераций в Soft_part_I
Results = struct('Begs',[],'Fins',[],'Types',[],'Centers',[], 'Widths',[],'Num',0,'n_iter',0);
Results = repmat(Results, Nstep, 1);

% Расчёт (или установка вручную) порога

% Показываем на графике несколько первых спектров
% figure(100);
% plot(abs( Y(:,1)),'b'); grid on; hold on
% plot(abs( Y(:,ceil(size( Y,2)*rand))),'r'); grid on; hold on
% plot(abs( Y(:,ceil(size( Y,2)*rand))),'g'); grid on; hold on
% % Запрашиваем выставление уровня порога
% [~,trh] = ginput(1);
% % close(100);
% fprintf(1,'Current trh = %g\n', trh);
% %trh = graythresh(abs(Y(:,:)));
% fprintf(1,'Current trh = %g\n', trh);

trh = 200;

Centras = zeros(Nstep, 1);
%% ==== Плисовская + I часть проц. алг-ма (первичная обработка) ==========
for k = 0:Nstep-1 %
    num = k;
    
    % ---- FPGA-часть алгоритма -----------------------------------------
    % [f,Sf,Sf_fltd,trh,des,Begs,Fins] = FPGA_part(samples, Nwindow, trh, TrhManual);

     Sf = Y(:, k+1); % + snr*nn;
    
    
  %% << вычисление ширины спектра  
    
    % Сглаживаем спектр, ибо слишком много "шерсти" в записях
    if Nwindow == 8192
        Nfilt = 10;
        Sf_fltd = conv(Sf,ones(Nfilt,1)/Nfilt,'same');
    else
        Nfilt = log2(Nwindow)-4; % Похоже, такой порядок оптимален
        Sf_fltd = conv(Sf,ones(Nfilt,1)/Nfilt,'same');
    end

 %    I = Sf_fltd > trh;
    I = Sf      > trh;
    
    des = zeros(size(Sf_fltd)); des(I) = 1.0;
    
    % Группировка превышений по критерию ширины
    % Если расстояние между превышениями меньше ширины меньшего из
    % превышений, то превышения сливаются
    dI  = diff(des);
    I10 = find(dI > 0); % Ищем все "фронты" (начала отрезков с превышен.)
    I01 = find(dI < 0); % Ищем все "спады" (концы отрезков с превыш.)

    if (isempty(I10) && isempty(I01))
            % Ничего не обнаружено, выдаём пустые границы
        Begs = [];
        Fins = [];
%         break
    end
    
    if(numel(I10) > 0)
    % Защита от ситуации, когда есть один фронт, но нет спада и наоборот
    if (length(I10) == length(I01)+1)
        % В конце спектра есть фронт, но нет спада - добавляем
        I01 = [I01; Nwindow/2];
    elseif (length(I10)+1 == length(I01))
        % В начале спектра есть спад, но нет фронта
        I10 = [1; I10];
    elseif (abs(length(I10)-length(I01)) > 1)
        % Такого вообще-то не должно быть, но на всякий случай
        error('Alghoritm error: diff of lengthes is more than 2!');
    else 
        % Длины векторов одинаковы
        if (I10(1) > I01(1))
            if (I10(end) > I01(end))
                % Произошёл "сдвиг" (началось спадом, закончилось фронтом)
                I10 = [1; I10];
                I01 = [I01; Nwindow/2];
            else
                % А это уже что-то странное!
                error('Strange situation!');
            end
        end
    end
    
    end % numel

    if numel(I10) > 0 % Вообще что-нибудь обнаружено?
        % Все превышения целиком находятся в текущем отрезке
        % (первый - фронт, последний - спад?)
        if (I10(1) < I01(1)) && (I10(end) < I01(end))
            Begs = I10; % Начала
            Fins = I01; % Концы
        else
            % Либо началось со спада, либо закончилось фронтом,
            % такого тоже не должно быть, но на всякий случай
            f0 = [1:length(Sf)];
            figure; % Технологическое окно
            % Выводим исходный спектр
            plot(f0,Sf,'g'); grid on; hold on
            % Выводим сглаженный спектр
            plot(f0,Sf_fltd,'k'); hold on
            % Показываем уровень порога
            plot([min(f0); max(f0)],[trh; trh],'r'); hold on
            disp(I10);
            disp(I01);
            error('Special situation!');
        end
    else
        % Ничего не обнаружено, выдаём пустые границы
        Begs = [];
        Fins = [];
    end
    
    
    %% ---- I часть процессорного алгоритма ------------------------------
    % Объединение превышений и предварительная оценка ширины спектра объединённых превышений
    Results(k+isMlab) = Soft_part_I(Sf,Begs,Fins,Width);
    if numel(Results(k+isMlab).Centers) > 0
    Centras(k+isMlab) = Results(k+isMlab).Centers(1);
    else
    Centras(k+isMlab) = NaN;    
    end
     % ---- Выводим рез-ты работы первых двух частей ---------------------
    figure(999);
    for tmp = 1:Results(k+isMlab).Num
        center = round(Results(k+isMlab).Centers(tmp));
        wdth   = Results(k+isMlab).Widths(tmp);
        plot([num+1 num+1],[Results(k+isMlab).Begs(tmp) Results(k+isMlab).Fins(tmp)],'g','LineWidth',2);
    end
    
end



% Centrs(:) = Centras(9:31);
% % Centrs(:) = Centras(2:23);
% x = 1:length(Centrs);
% 
% 
% 
% c = polyfit(x,Centrs,1);
% disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% 
% y_est = polyval(c,x);
% 
% figure
% plot(x, Centrs, '.-b', x,y_est,'.-r')
% legend('center freq', 'fitted line')
% xlabel('Nstep')
% ylabel('Bins')
% grid on
% 
% disp(['Error is ' num2str(sum((y_est-Centrs).^2)/length(x))]);

% return
% ---- Проверка на заниженный порог --------------------------------------
tmp = mean([Results(:).Num]); % Среднее число отдельных отметок
if (tmp > MaxMeanPlts)
    xlabel('Too many plots detected, raise the threshold','Color','r');
    warning('Too many plots detected, raise the threshold');
    return
end

%% ==== Поиск и построение трасс  ========================================
deltas = [delta1 delta2]; % "Зоны" захвата и сопровождения (см. далее)

% Nplots  = 0; % Общее кол-во превышений

% % % Для всех отметок создаётся дополнительные поля isInTrace и SerNumber, 
% % % при этом: 
% % %    -1 означает, что отметка бесперспективна (исходно для всех);
% % %     0 означает, что никуда не входит;
% % %     n > 0 означает, что отметка входит в трассу с номером n, 
% % % а также все отметки нумеруются подряд сверху вниз и слева направа
% % % (в текущей версии нумерация использ. только для удобства)
for i = 1:Nstep 
    Results(i).isInTrace = -1*ones(size(Results(i).Centers)); 
    %Results(i).SerNumber = zeros(size(Results(i).Centers)); 
end

% % Нумеруем подряд все отметки
% for i = 1:Nstep 
%     %Results(i).SerNumber = [Nplots+1 : Nplots+Results(i).Num]; 
%     Nplots = Nplots+Results(i).Num;
% end

%  ===== Предварительный отсев бесперспективных отметок ===================
for i = 1:Nstep-2
    for k = 1:Results(i).Num % Перебор всех отметок (превышений)
        %{
    Есть возможность продолжить текущее превышение в соотвествии
    с выбранным критерием?
    Текущий критерий - если в пределах длины текущей отметки
     +/- deltas[1] от длины есть следующая отметка, или следующей нет, но
     +/- deltas[2] от длины есть отметка через одну
        %}
        Begn = Results(i).Begs(k); % Начало текущей отметки
        Finh = Results(i).Fins(k); % Конец текущей отметки
        Widh = Results(i).Widths(k); % Протяжённость текущей отметки
        
        % Сначала исследуем на один шаг вперёд
        for k1 = 1:Results(i+1).Num
            Begn1 = Results(i+1).Begs(k1); % Начало текущей отметки
            Finh1 = Results(i+1).Fins(k1); % Конец текущей отметки
            Widh1 = Finh1-Begn1; % Протяжённость текущей отметки
            % Отметки перекрываются с учётом допуска (deltas(1))?
            % Прим.: допуск, вообще говоря, должен считаться по размеру БПФ
            % и текущему ОСШ
            if (Begn*(1-deltas(1)) > Finh1) || (Finh*(1+deltas(1)) < Begn1)
                % Текущая отметка с номером k1 не проходит - переходим к следующей
                continue
            else
                % Есть чем продолжить - отметка переходит в нулевой класс
                Results(i).isInTrace(k) = 0;  break
            end
        end
        
        if (Results(i).isInTrace(k) == 0), continue, end % К следующему k или i
        
        % Если не найдено ничего среди соседей - смотрим на два шага вперёд
        for k1 = 1:Results(i+2).Num
            Begn1 = Results(i+2).Begs(k1); % Начало текущей отметки
            Finh1 = Results(i+2).Fins(k1); % Конец текущей отметки
            Widh1 = Finh1-Begn1; % Протяжённость текущей отметки
            % Отметки перекрываются с учётом допуска (deltas(2))?
            if (Begn*(1-deltas(2)) > Finh1) || (Finh*(1+deltas(2)) < Begn1)
                % Текущая отметка с номером k1 не проходит - переходим к следующей
                continue
            else
                % Есть чем продолжить - отметка переходит в нулевой класс
                Results(i).isInTrace(k) = 0;  break
            end
        end
        
        % Если и после второго шага не найдено отметки, которая
        % перекрывается с текущей, то Results(i).isInTrace(k) остаётся
        % равным -1 и отметка выбрасывается из дальнейшего рассмотрения
        
        % Надо бы написать более универсальную процедуру, в которой кол-во
        % "дырок" было бы параметром
    end
end

% ---- Только для Матлаба - выводим результат ----------------------------
figure(999);
for k = 1:Nstep-2
    for tmp = 1:Results(k).Num
        if Results(k).isInTrace(tmp) >= 0
            center = round(Results(k).Centers(tmp));
            wdth   = Results(k).Widths(tmp);
            % Выводим "выжившие" отметки красным цветом
            plot([k k], [Results(k).Begs(tmp) Results(k).Fins(tmp)],'r','LineWidth',2);
        end
    end
end
% ------------------------------------------------------------------------

% return

% ===== Завязываем трассы ====================================
% Тактика такова: берём очередную отметку, у которой
% Results(i).isInTrace(k) = 0 и пытаемся от неё продлить трассу так, что бы
% этот "зачаток" будущей трассы (временная трасса) имел длину не менее Ncptrd.
% При проводке трассы мы, при наличии альтернатив, выбираем в кач-ве след.
% отметки ту, центр которой максимально близок к центру текущей отметки (весьма дискутируемо!)
% Если упираемся в отметку, которая уже входит в какую-либо трассу, то (в
% текущей версии) проводка трассы прекращается и все ранее объёдинённые во
% временную трассу отметки получают тип -1 дабы не плодить тучу ветвящихся трасс
% (тоже дискутируемо)
% Если это удаётся - то новая трасса получает подтверждение и очередной номер и мы
% продолжаем её тянуть пока не нарвёмся на отметку с типом -1. 
% Подтверждённую трассу продолжаем проводить даже если он попадает на
% отметки, которые уже входят в ранее проведённую трассу, при этом для
% таких отметок isInTrace++. 
% По идее, при такой стратегии проводки, автоматически произойдёт защита от
% ШП-помех (трасса пройдёт просто через них), но, возможно, для "узкоп."
% трасс имеет смысл ввести отдельную стартегию проводки с калмановской
% фильтрацией (хотя, если "излом" трассы будет закрыт ШП-импульсом, это
% может быть фатально - сопровождение будет сорвано)
% Если не удаётся набрать Ncptrd отметок во временной трассе, то все
% отметки из временной трассы получают -1 (хотя, возможно, стоило бы
% откатиться назад к последней отметке, имевшей альтернативное продолжение
% проводки, но пока делать этого не будем)

% Резервируем для временной трассы
traces = struct('Res_i',0,'Num_in_Res',0,'Center',0,'Clutter',0);
traces = repmat(traces, Nstep, MaxTracks); % Колонка - одна трасса

Ntraces = 0; % Кол-во обнаруженных "треков"

flag_too_many_traces = 0; % Для завершения поиска при перегрузке трассами

% Перебор отметок, с которых можно начать врЕменную трассу
for i = 1:Nstep-2
    for k = 1:Results(i).Num
        
        if (Results(i).isInTrace(k) == 0)
            
            % "Обнуляем" временную трассу
            temp_trace  = struct('Res_i',0,'Num_in_Res',0,'Center',0,'Clutter',0);
            temp_trace  = repmat(temp_trace, Nstep, 1);
            temp_length = 1;
            isClutter   = 0;

            % Вектор для расчёта медианной ширины отметок в трассе
            %temp_lengs = zeros(Nstep, 1); 
            
            % Первая отметка - текущая
            temp_trace(1).Res_i      = i;
            temp_trace(1).Num_in_Res = k;
            temp_trace(1).Center     = Results(i).Centers(k);
            %temp_lengs(1)            = Results(i).Widths(k);
            
            i2 = i; k2 = k; 
            
            % Теперь пытаемся продолжить врЕменную трассу
            % перебор всех неотброшенных и не входящих в другие трассы отметок
            while(1)
                % Находим ближайшую к последней найденной
                if ((temp_length > 1) && (isClutter == 1))
                    isClutter_prev = isClutter;
                else
                    isClutter_prev = 0;
                end
                i_prev = i2; k_prev = k2;
                [tmp,k2,isClutter] = GetNextClosest(Results(i_prev),k_prev,Results(i2+1:i2+2));
                i2 = i2+tmp;
                if (tmp > 0) % Что-нибудь найдено?
                    %fprintf(1,'Temp trace #i %g, k = %g, k2 = %g, Results(i).isInTrace(k) = %g\n', i, k, k2, Results(i).isInTrace(k));
                    if (Results(i2).isInTrace(k2) > 0)
                        % Упёрлись в отметку, которая уже куда-то входит - помечаем все отметки из
                        % временной трассы -1 и уходим на поиск следующих трасс
                        % Здесь, кстати, дискутируемо!
                        Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                        break
                        
                    elseif (temp_length >= Ncptrd)
                        % Трасса подтверждена - проверяем не превышено ли допустимое кол-во трасс
                        if (Ntraces >= MaxTracks)
                            warning('Too many traces have been detected!'); %#ok<WNTAG>
                            % Слишком много трасс - поиск закончен
                            flag_too_many_traces = 1; 
                            break
                        else
                            % Иначе - создаём новую трассу
                            Ntraces = Ntraces+1; % Новая трасса
                            traces(:,Ntraces) = temp_trace;
                            Results = SetWholeTraceTo(Results,temp_trace,temp_length,Ntraces);
                            
                            % И продолжаем трассу до куда возможно
                            while(1)
                                % Находим ближайшую к последней найденной
                                isClutter_prev = isClutter; % Задержка признака на один шаг
                                if ((isClutter == 1))
                                    % Последняя отметка пришла с признаком "помеха" - нужно её перепрыгнуть
                                    i_prev = i2; % traces(tmp-1,Ntraces).Res_i; % Просто i2-1 нельзя - можем сорвать слежение
                                    k_prev = k2; %traces(tmp,Ntraces).Num_in_Res; % Та же фигня
                                    %i2 = i2+1; % Перепрыгиваем через импульсную помеху
                                else
                                    i_prev = i2; k_prev = k2;
                                end
                                [tmp,k2,isClutter] = GetNextClosest(Results(i_prev),k_prev,Results(i2+1:i2+2));
                                i2 = i2+tmp; % На сколько шагов смещаемся (1 или 2, 2 или 3, если предыдущий был помехой)
                                if (tmp > 0)
                                    % Ищем первое свободное (нулевое) место
                                    tmp = find(~[traces(:,Ntraces).Res_i],1); % Удлиняем трассу
                                    traces(tmp,Ntraces).Res_i       = i2;
                                    traces(tmp,Ntraces).Num_in_Res  = k2;
                                    traces(tmp,Ntraces).Center      = Results(i2).Centers(k2); % Здесь нужно подправить, если isClutter (см. далее)
                                    traces(tmp,Ntraces).Clutter     = isClutter; % Признак "помеха"
                                    Results(i2).isInTrace(k2)       = Ntraces;
                                    % Подправляем предыдущую отметку, если она пришла с признаком "помеха", иначе трасса может уйти в сторону 
                                    if (isClutter),
                                        x = [traces(tmp-Nlin : tmp-1, Ntraces).Res_i];
                                        Y = [traces(tmp-Nlin : tmp-1, Ntraces).Center];
                                        xi = traces(tmp, Ntraces).Res_i;
                                        % вместо фактического центра подставляем линейную экстраполяцию по Nlin предыдущим отметкам
                                        traces(tmp,Ntraces).Center = interp1(x, Y, xi,'linear','extrap');
                                    end
                                    
                                    % ---- Только для отладки --------------------------------------
                                    fprintf(1,'================================\n');
                                    fprintf(1,'x = %g\n', traces(tmp,Ntraces).Res_i);
                                    fprintf(1,'y = %g\n', traces(tmp,Ntraces).Num_in_Res);
                                    fprintf(1,'Center = %g\n', traces(tmp,Ntraces).Center);
                                    fprintf(1,'isClutter %g\n', traces(tmp,Ntraces).Clutter);
                                    fprintf(1,'Trace No.%g\n', Ntraces);
                                    % ------------------------------------------------------------------------
                                    
                                    % Временную трассу тоже продолжаем тащить
                                    temp_length                        = temp_length+1; % Увеличили длину трассы
                                    temp_trace(temp_length).Res_i      = i2;
                                    temp_trace(temp_length).Num_in_Res = k2;
                                    temp_trace(temp_length).Center     = traces(tmp,Ntraces).Center; % и здесь!!!
                                    temp_trace(temp_length).Clutter    = isClutter;
                                    %temp_lengs(temp_length)            = Results(i2).Widths(k2);

                                    if (i2 >= length(Results)-3), % Что бы не вылететь за пределы спектрограммы 
                                        ShowMeTheTrace(Results,traces,Ntraces);
                                        break; 
                                    end
                                else
                                    % Нечем продолжать - рисуем трассу и переходим к поиску новых временных трасс
                                    %Results = SetWholeTraceTo(Results,traces(:,Ntraces),tmp,Ntraces);
                                    ShowMeTheTrace(Results,traces,Ntraces);
                                    %pause;
                                    break % Трасса с номером Ntraces завершена
                                end
                            end
                            
                            % Уходим на поиск следующей трассы
                            break
                        end
                        
                    else
                        % Нормальное продолжение временной трассы
                        temp_length                        = temp_length+1;
                        temp_trace(temp_length).Res_i      = i2;
                        temp_trace(temp_length).Num_in_Res = k2;
                        temp_trace(temp_length).Clutter    = isClutter;
                        %temp_lengs(temp_length)              = Results(i2).Widths(k2);
                        temp_trace(temp_length).Center     = Results(i2).Centers(k2); % И тут тоже!
                        
                        if (i2 >= length(Results)-2),
                            % ВрЕменной трассе не повезло - спектрограмма
                            % кончилась
                            Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                            break;
                        end
                        
                    end
                else
                    % Нет никакого продолжения - помечаем все отметки из
                    % временной трассы -1 и уходим на поиск следующих
                    % трасс
                    Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                    %flag_to_exit = 1; 
                    break
                end
                
            end % of while
            
            if (flag_too_many_traces), break, end
        end
        % Конец цикла проводки врЕменной трассы
        
        
    end
    if (flag_too_many_traces), break, end
end

% ---- Только для Матлаба - выводим результат ----------------------------
% figure(999);
% for i = 1:Ntraces,
%     ShowMeTheTrace(Results, traces, i);%
% end
% figure(999); xlabel('Classified: RAMERA type radar','FontSize', 16);
% ------------------------------------------------------------------------

% ==== Оцениваем заданные параметры треков ===============================
% Оцениваем: среднюю частоту, максимальное изменение частоты, наличие
% "изломов", линейный и квадратичный члены полиномиальной апроксиммации
% трека (насколько он отличается от ЛЧМ), СКО частоты, отношение СКО
% частоты трека к длине трека (избавление от треков-"туманностей"),
% длину трека, etc.

% Во-первых, проверяем число обнаруженных треков, если их слишком много,
% то явно был неверный выбор порога, либо слишком малое ОСШ для
% распознования. Такое эпизод - отбраковывается.
if (Ntraces > MaxTraces)
    fprintf(1,'===== Final Result ====================================\n');
    fprintf(1,'Too many traces: unclassified.\n');
    return
end

% Сюда пишем оценЁнные параметры треков
TracksInfo = struct('Fmean',0,'dFmax',0,'Length',0,'Curvature',0,'MSE',0,'isBreak',0);
TracksInfo = repmat(TracksInfo, Ntraces);
Lengths    = zeros(Ntraces, 1); % Длины трасс

for n_trk = 1:Ntraces
    the_last = find(~[traces(:,n_trk).Res_i],1)-1; % Длина текущей трассы
    Lengths(n_trk) = the_last;
    
    x        = [traces(1:the_last,n_trk).Res_i];
    tmp      = [traces(1:the_last,n_trk).Num_in_Res];
    wtdhs    = zeros(size(x));
    
    for i = 1:length(x)
        wtdhs(i) = Results(x(i)).Types(tmp(i)); % Типы отметок по широкополосности
    end
    
    % Проверка на Сталкер - единственный пока в базе ШП-радар
    if (mean(wtdhs)/length(wtdhs) > 0.5),
        fprintf(1,'===== Final Result ====================================\n');
        fprintf(1,'Most probable radar: %s\n', 'Stalker 34G');
        fprintf(1,'With probability: %g%%\n', mean(wtdhs)/length(wtdhs));
        return % Выходим - дальше исследовать незачем, других ШП в базе пока нет
    end
    
    % Собираем трек целиком
    y = [traces(1:the_last,n_trk).Center]; % Центры отметок
    
    % Замеряем основные параметры
    TracksInfo(n_trk).Fmean  = mean(y);
    TracksInfo(n_trk).dFmax  = max(y)-min(y);
    TracksInfo(n_trk).Length = x(end)-x(1);
    p = polyfit(x,y,2); TracksInfo(n_trk).Curvature = p(1)/p(2);
    p = polyfit(x,y,10); Y = polyval(p,x); err = abs(Y-y);
    TracksInfo(n_trk).MSE = std(err);
    
    % Проверка на наличие "изломов" треков
    trh_brk = mean(err)+3.5*std(err); % Порог обнаружения изломов
    
    if ~isempty(find(err > trh_brk,1))
        % Трек с изломом - яркий классификационный признак
        TracksInfo(n_trk).isBreak = 1; % Обнаружен излом трека
        [~,i_max] = max(err); % Засекаем положение излома трека
    end
    
end

% ---- Только для Матлаба - выводим результат распознавания --------------
Types = {'Non-radar signal','Iskra','MultaNova','Multa-CD',...
    'MultaNova6F','Ramera-222','Stalker 34G'};
fprintf(1,'===== Final Result ====================================\n');

fprintf(1,'Traces detected: %g\n', Ntraces);
%fprintf(1,'Possible to comb.: 0'); 
for n_trk = 1:Ntraces
    fprintf(1,'---- Track No.%g -----------------\n',n_trk);
    fprintf(1,'F_{mean}  = %g\n', TracksInfo(n_trk).Fmean); 
    fprintf(1,'dF_{max}  = %g\n', TracksInfo(n_trk).dFmax); 
    fprintf(1,'Curvature = %g\n', TracksInfo(n_trk).Curvature); 
    fprintf(1,'MSE       = %g\n', TracksInfo(n_trk).MSE); 
    if (TracksInfo(n_trk).isBreak)
        fprintf(1,'Fracture(s): YES\n');
    else
        fprintf(1,'Fracture(s): NO\n');
    end
end

% ---- Собственно распознавание и вывод результатов ----------------------
if (Ntraces == 1)
    if (TracksInfo(1).dFmax < 10)
        % Что-то узкополосное (возможно "Мультанова 6F", но увы)
        Type_of_radar = 'Unclassified';
        Prob_of_recogn = 0.0;
    else
        if (TracksInfo(1).isBreak == 1)
            Type_of_radar = 'Ramera-222';
            Prob_of_recogn = TracksInfo(n_trk).dFmax/(Nwindow/8);
        else
            if (TracksInfo(1).dFmax > 200)
                Type_of_radar = 'Ramera-222';
                Prob_of_recogn = 1-1/(TracksInfo(1).dFmax/25);
            else
                Type_of_radar = 'Unclassified';
                Prob_of_recogn = 0.0;
            end
        end
    end
else
    % Если треков больше, то главный вопрос насколько они перекрываются
    if (IntersectEnough(traces(:,1:Ntraces),Lengths))
        % Если они достаточно длинные и хорошо перекрываются - то это,
        % скорее всего "Мульта-CD"
        Type_of_radar = 'Multa CD';
        tmp = [TracksInfo(:).Length];
        Prob_of_recogn = min(tmp(1:2))/max(tmp(1:2));
    else
        if ((Ntraces == 2) || (Ntraces == 3))
            % Если треки не перекрываются и треков два или три, то это,
            % скорее всего "Рамера", но зависит от того как далеко треки
            % друг от друга
            Type_of_radar = 'Iskra';
            tmp = [TracksInfo(:).Length];
            Prob_of_recogn = (traces(1,2).Res_i-traces(Lengths(1),1).Res_i)/Nstep;
        else
            % Что-то не понятное
            Type_of_radar = 'Unclassified';
            Prob_of_recogn = 0.0;
        end
            
    end
end

fprintf(1,'Most probable radar: %s\n', Type_of_radar);
fprintf(1,'With probability: %g%%\n', Prob_of_recogn*100);

figure(999);
xlabel(sprintf('Most probable radar: %s with P_D=%g', Type_of_radar, Prob_of_recogn*100),'Color','r');
% ------------------------------------------------------------------------



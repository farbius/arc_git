function [ isIntersect ] = IntersectEnough(traces,lengths)
% ƒе-юре функци€ определ€ет имеетс€ ли в списке обнаруженных трасс traces с
% длинами, задаваемыми вектором lengths,
% хот€ бы пара досаточно длинных трасс и, если имеетс€, то перекрываютс€ ли
% по оси x (номер окна) эти трассы в достаточной степени.
% ƒе-факто, в данном релизе программы, определ€ет сигнал "ћульты-CD", т.к.
% в запис€х нет друго радара с несколькими узкополосными треками (пока)

global Nstep
% ---- ћестные параметры алгоритма ---------------------------------------
LimLength   = 0.9; % ћинимальна€ потребна€ длина трека к длине окна
LimIntersec = 0.5; % ћинимальное перекрытие треков в дол€х длины

Ntraces = size(traces,1);

Ilm = find(lengths/Nstep >= LimLength);
if length(Ilm) < 2, % ≈сть хот€ бы два достаточно длинных трека?
    isIntersect = 0;
else
    Ilm
    lengths(Ilm)
    size(traces)
    x1 = traces(1,Ilm(1)).Res_i;
    x2 = traces(lengths(Ilm(1)),Ilm(1)).Res_i;
    x3 = traces(1,Ilm(2)).Res_i;
    x4 = traces(lengths(Ilm(2)),Ilm(2)).Res_i;
    tmp = min(x4,x2) - max(x3,x1);
    
    if (tmp/Nstep > LimIntersec),
        isIntersect = 1;
    else
        isIntersect = 1;
    end
end


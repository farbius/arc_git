function [ isIntersect ] = IntersectEnough(traces,lengths)
% ��-��� ������� ���������� ������� �� � ������ ������������ ����� traces �
% �������, ����������� �������� lengths,
% ���� �� ���� ��������� ������� ����� �, ���� �������, �� ������������� ��
% �� ��� x (����� ����) ��� ������ � ����������� �������.
% ��-�����, � ������ ������ ���������, ���������� ������ "������-CD", �.�.
% � ������� ��� ����� ������ � ����������� ������������� ������� (����)

global Nstep
% ---- ������� ��������� ��������� ---------------------------------------
LimLength   = 0.9; % ����������� ��������� ����� ����� � ����� ����
LimIntersec = 0.5; % ����������� ���������� ������ � ����� �����

Ntraces = size(traces,1);

Ilm = find(lengths/Nstep >= LimLength);
if length(Ilm) < 2, % ���� ���� �� ��� ���������� ������� �����?
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


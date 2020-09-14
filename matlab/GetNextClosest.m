function [i2,n2,isClutter] = GetNextClosest(Results1,k,Results2)
% ������� �������������� ����� ���������� ������ ��� ������� ������� 
% (Results1), � ��������� ����(Results2(1)), ��� ����� ���� (Results2(2))
% ���� �� ��������� �������� � ������� deltas
% ���� ��������������, ��� ����� ��������� ������������ ������� �
% clutterTrh ��� ��� ����� ����������� ����� ������� �������, �� ���������
% ������� ����������� ���������� ������ (isClutter = 1).

%
% ����������:
%   i2  - �������� ��������� ������� (0 - �� ������� �����������)
%   n2  - ����� ���������� � ������� (�� 1 �� Results2(i2).Num)

global deltas clutterTrh

i2 = 0; n2 = 0; % ���������� ������ �� �������

isClutter = 0; % ����� ���� ���� �� �����

found = []; % ���� ���������� ����������

Begn = Results1.Begs(k); % ������ ������� �������
Finh = Results1.Fins(k); % ����� ������� �������
Wdth = Results1.Widths(k); % ������ ������� �������

% ������� ��������� �� ���� ��� �����
for k1 = 1:Results2(1).Num,
    if (Results2(1).isInTrace(k1) == 0),
        Begn1 = Results2(1).Begs(k1); % ������ ������� �������
        Finh1 = Results2(1).Fins(k1); % ����� ������� �������
        % ������� ������������� � ������ ������� (deltas(1))?
        if (Begn*(1-deltas(1)) > Finh1) || (Finh*(1+deltas(1)) < Begn1),
            % ������� ������� � ������� k1 �� �������� - ��������� � ���������
        else
            % ���� ��� ���������� - ������� ������� � found
            found = [found; k1];
        end
    end
end

if ~isempty(found),
    i2 = 1; 
    if length(found) == 1,
        n2 = found; % ���������� ����� ������ ����� �������� - �� ������.
    else
        % ����� ���������� ����������� ���������, ����� �� �������:
        % �������� ��������� �� ����� � ��������� ����� ������
        %[~,n2] = min(abs(Results2(1).Widths - med_len);
        % �������� ��������� �� ����� � ����� ����������
        %[~,n2] = min(abs(Results2(1).Widths - Results1.Widths(k))); -- ��
        %������ ��������!
        [~,n2] = min(abs(Results2(1).Centers - Results1.Centers(k)));
    end
    if (Results2(1).Widths(n2)/Wdth >= clutterTrh),
        isClutter = 1; % ���� ����������, ��� ���������� �������
    end
    return % ����� ��������
end

% ���� �� ������� ������ ����� ������� - ������� �� ��� ���� �����
for k1 = 1:Results2(2).Num,
    if (Results2(2).isInTrace(k1) == 0),
        Begn1 = Results2(2).Begs(k1); % ������ ������� �������
        Finh1 = Results2(2).Fins(k1); % ����� ������� �������
        % ������� ������������� � ������ ������� (deltas(2))?
        if (Begn*(1-deltas(2)) > Finh1) || (Finh*(1+deltas(2)) < Begn1),
            % ������� ������� � ������� k1 �� �������� - ��������� � ���������
        else
            % ���� ��� ���������� - ������� ������� � found
            found = [found; k1];
        end
    end
end

if ~isempty(found),
    i2 = 2;
    if length(found) == 1,
        n2 = found;
    else
        % ����� ���������� ����������� ���������,����� -
        % �������� ��������� �� ������������ (Center)
        % �� ����, ������, ������ �������� ������������ "������"
        [~,n2] = min(abs(Results2(2).Centers - Results1.Centers(k)));
    end
    if (Results2(2).Widths(n2)/Wdth >= clutterTrh),
        isClutter = 1; % ���� ����������, ��� ���������� �������
    end
end

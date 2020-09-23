%�����������          : ?
%���� ����.           : 08.09.2015
%���� ������          : 22.04.2016
%������               : 1.2	    
%������������         : �������, ����������� ����� ���������,
%                       ��������������� ��� �������� � ���������
%�����������          : ����� I ������ ���������� ������ ��� �� ��������� 
%                       ������ ������ ������ �� FPGA � ������ ���������� 
%                       �� ����� 41���. 
function [Results] = Soft_part_I(Sf,Begs,Fins,Width)

isMlab = 1;
Nlength = length(Sf);
% ===== ��������� FPGA-��������� =========================================
% ���������� ����������
if (any(Begs) == 0)  % ���������� ����������?
    % ���� ���, �� ��������� ������ �������� ��� ���������� ����
    Results.Begs    = [];
    Results.Fins    = [];
    Results.Types   = [];
    Results.Centers = [];
    Results.Widths  = [];
    Results.Num     = 0;
    Results.n_iter  = 0;
    return % ������� �� �������
end

% ==== ��������������� ��������� ���������� (����������� ����������) =
n_iter = 0; % ��������� ������������!
while 1
    Lens = Fins-Begs; % ����� �������� ����������
    if (length(Lens) == 1)
        break; % ���� ��������� ������������ ��������, �� ������� �� while
    end
    Bets   = zeros(length(Lens)-1,1); % ������ ���������� ����� ��������� ����������
    IsJoin = zeros(length(Lens)-1,1); % ���� ���������� �����������
    for i = 1:length(Lens)-1
        %         Bets(i) = Begs(i+1)-Fins(i); % ����� �����������
        %         % �������� �����������:
        %         % ���� ���������� ����� ��������� ��������� ������ ������
        %         % ������ �� ��������, �� ���� �������� ������������ � ����
        %         if ((Begs(i+1) > 0) && (Fins(i) > 0))
        %             if ((Lens(i) >= Bets(i))||(Bets(i) <= Lens(i+1))),
        %                 IsJoin(i) = 1; % �������� ���������� �� ��������
        %             end
        %         end
        
        Spare = 5.0; % ��� �� �� "�������" �� ����� ������� ����������
        Bets(i) = Spare*(Begs(i+1)-Fins(i)); % ����� �����������
        % ���� ���������� ����� ��������� ��������� ������ ������
        % ����� �������� � �������, �� ���� �������� ������������ � ����
        if ((Begs(i+1) > 0) && (Fins(i) > 0))
            if ((Lens(i) >= Bets(i)) || (Bets(i) <= Lens(i+1)))
                IsJoin(i) = 1; % �������� ���������� �� ��������
            end
        end
    end
    
    if (sum(IsJoin) == 0)
        break; % ���� ������ ������ ����������, �� ������� �� while
    end
    
    % ���������� ����������� ��������
    I = find(IsJoin == 1);
    Begs(I+1) = []; % ������� ������ � �������� I+1
    Fins(I)   = []; % ������� ����� � �������� I
    
    n_iter = n_iter+1;
    if (n_iter > 50)
        % ������� ����� ����������, ������ - ���, �� ���� ������ �������
        warning('The number of iterations exceeds 20!');
        break;
    end
end % ������ �� ��������� ��������

% ������ ������� ������� � ������ ���������� ����������
Centers = zeros(size(Begs));
Widths  = zeros(size(Begs));
for i = 1:length(Begs)
    x = Begs(i):Fins(i);
    Centers(i) = sum(x'.*Sf(x).^2)/sum(Sf(x).^2);
    %Widths(i)  = sum((x'-Centers(i)).^2.*Sf(x).^2)/sum(Sf(x).^2);
    Widths(i)  = length(x);
    if (Widths(i) > Nlength)
        Widths(i) = Nlength; % ��������� ������� ��� ��� ������ �������!
    end
end

% ������������� �� ����������������
% I = ((Fins-Begs) > Width);
I = (Widths > Width);

% ��������� �������� ��� �������� ����
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

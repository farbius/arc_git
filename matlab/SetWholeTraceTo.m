function Results = SetWholeTraceTo(Results,temp_trace,temp_length,num)
% ������� ������������� ��� �������, �������� �� ��������� ������
% temp_trace � �������� num �������� �������� isInTrace
for i = 1:temp_length
    Results(temp_trace(i).Res_i).isInTrace(temp_trace(i).Num_in_Res) = num;
end

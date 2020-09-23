function Results = SetWholeTraceTo(Results,temp_trace,temp_length,num)
% ‘ункци€ предназанчена дл€ отметок, вход€щих во вр≈менную трассу
% temp_trace в заданный num значение признака isInTrace
for i = 1:temp_length
    Results(temp_trace(i).Res_i).isInTrace(temp_trace(i).Num_in_Res) = num;
end

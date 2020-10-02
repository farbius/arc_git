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
% ��������� ���������� 
% ========================================================================
Fs          = 800e6; %| ������� ���������
Nwindow     = 8192;  %| �������� ���� ���
Width       = 40;    %| ����� �������� ������� � "����������������" (� ������ ���)
TrhManual   = 1;     %| ������ ��������� ������ (1 - ��, 0 - �������)
MaxMeanPlts = 7.0;   %| ������������ ������� ���-�� ��������� ����������
delta1      = 0.25;  %| ������������ ����������� �������� ������� �� 1 ���
delta2      = 0.5;   %| ������������ ����������� �������� ������� �� 2 ����
Ncptrd      = 30;    %| ����������� ���-�� ������� ��� ������� ������
MaxTracks   = 30;    %| ������������ ���-�� ����� (�� ��������� ����������)
MaxTraces   = 5;     %| ������������ ���-�� ����� ��� �������������
clutterTrh  = 3.0;   %| ��. ��������� GetNextClosest
Nlin = 3; %| ���-�� �������, �� ������� �������������� ��������� �������
% ========================================================================
% ����������: ������ ������, delta1 � delta2 ������ �����������
% ������������ �� �������� �������� ��������� � �������� ���. ������, ���
% ��������� ���������� ��������� ������� � ������� ��-������� ���
% ����������� "�������" ������������ ������� (��� ������ ������ ���������
% Soft_part_I, �� �� �������������� ���������)
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

% ���� �� ��������������� ����� - ����������� ���� ��������
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


% ��������� ��� ������� ���� � ������������ ����������,
% ����������� � �������������
% ����: ������ �������� � �����������, ����� -//-, ���, ������ ��������,
% ����� ���-�� ��������, ���-�� �������� � Soft_part_I
Results = struct('Begs',[],'Fins',[],'Types',[],'Centers',[], 'Widths',[],'Num',0,'n_iter',0);
Results = repmat(Results, Nstep, 1);

% ������ (��� ��������� �������) ������

% ���������� �� ������� ��������� ������ ��������
% figure(100);
% plot(abs( Y(:,1)),'b'); grid on; hold on
% plot(abs( Y(:,ceil(size( Y,2)*rand))),'r'); grid on; hold on
% plot(abs( Y(:,ceil(size( Y,2)*rand))),'g'); grid on; hold on
% % ����������� ����������� ������ ������
% [~,trh] = ginput(1);
% % close(100);
% fprintf(1,'Current trh = %g\n', trh);
% %trh = graythresh(abs(Y(:,:)));
% fprintf(1,'Current trh = %g\n', trh);

trh = 200;

Centras = zeros(Nstep, 1);
%% ==== ���������� + I ����� ����. ���-�� (��������� ���������) ==========
for k = 0:Nstep-1 %
    num = k;
    
    % ---- FPGA-����� ��������� -----------------------------------------
    % [f,Sf,Sf_fltd,trh,des,Begs,Fins] = FPGA_part(samples, Nwindow, trh, TrhManual);

     Sf = Y(:, k+1); % + snr*nn;
    
    
  %% << ���������� ������ �������  
    
    % ���������� ������, ��� ������� ����� "������" � �������
    if Nwindow == 8192
        Nfilt = 10;
        Sf_fltd = conv(Sf,ones(Nfilt,1)/Nfilt,'same');
    else
        Nfilt = log2(Nwindow)-4; % ������, ����� ������� ���������
        Sf_fltd = conv(Sf,ones(Nfilt,1)/Nfilt,'same');
    end

 %    I = Sf_fltd > trh;
    I = Sf      > trh;
    
    des = zeros(size(Sf_fltd)); des(I) = 1.0;
    
    % ����������� ���������� �� �������� ������
    % ���� ���������� ����� ������������ ������ ������ �������� ��
    % ����������, �� ���������� ���������
    dI  = diff(des);
    I10 = find(dI > 0); % ���� ��� "������" (������ �������� � ��������.)
    I01 = find(dI < 0); % ���� ��� "�����" (����� �������� � ������.)

    if (isempty(I10) && isempty(I01))
            % ������ �� ����������, ����� ������ �������
        Begs = [];
        Fins = [];
%         break
    end
    
    if(numel(I10) > 0)
    % ������ �� ��������, ����� ���� ���� �����, �� ��� ����� � ��������
    if (length(I10) == length(I01)+1)
        % � ����� ������� ���� �����, �� ��� ����� - ���������
        I01 = [I01; Nwindow/2];
    elseif (length(I10)+1 == length(I01))
        % � ������ ������� ���� ����, �� ��� ������
        I10 = [1; I10];
    elseif (abs(length(I10)-length(I01)) > 1)
        % ������ ������-�� �� ������ ����, �� �� ������ ������
        error('Alghoritm error: diff of lengthes is more than 2!');
    else 
        % ����� �������� ���������
        if (I10(1) > I01(1))
            if (I10(end) > I01(end))
                % ��������� "�����" (�������� ������, ����������� �������)
                I10 = [1; I10];
                I01 = [I01; Nwindow/2];
            else
                % � ��� ��� ���-�� ��������!
                error('Strange situation!');
            end
        end
    end
    
    end % numel

    if numel(I10) > 0 % ������ ���-������ ����������?
        % ��� ���������� ������� ��������� � ������� �������
        % (������ - �����, ��������� - ����?)
        if (I10(1) < I01(1)) && (I10(end) < I01(end))
            Begs = I10; % ������
            Fins = I01; % �����
        else
            % ���� �������� �� �����, ���� ����������� �������,
            % ������ ���� �� ������ ����, �� �� ������ ������
            f0 = [1:length(Sf)];
            figure; % ��������������� ����
            % ������� �������� ������
            plot(f0,Sf,'g'); grid on; hold on
            % ������� ���������� ������
            plot(f0,Sf_fltd,'k'); hold on
            % ���������� ������� ������
            plot([min(f0); max(f0)],[trh; trh],'r'); hold on
            disp(I10);
            disp(I01);
            error('Special situation!');
        end
    else
        % ������ �� ����������, ����� ������ �������
        Begs = [];
        Fins = [];
    end
    
    
    %% ---- I ����� ������������� ��������� ------------------------------
    % ����������� ���������� � ��������������� ������ ������ ������� ����������� ����������
    Results(k+isMlab) = Soft_part_I(Sf,Begs,Fins,Width);
    if numel(Results(k+isMlab).Centers) > 0
    Centras(k+isMlab) = Results(k+isMlab).Centers(1);
    else
    Centras(k+isMlab) = NaN;    
    end
     % ---- ������� ���-�� ������ ������ ���� ������ ---------------------
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
% ---- �������� �� ���������� ����� --------------------------------------
tmp = mean([Results(:).Num]); % ������� ����� ��������� �������
if (tmp > MaxMeanPlts)
    xlabel('Too many plots detected, raise the threshold','Color','r');
    warning('Too many plots detected, raise the threshold');
    return
end

%% ==== ����� � ���������� �����  ========================================
deltas = [delta1 delta2]; % "����" ������� � ������������� (��. �����)

% Nplots  = 0; % ����� ���-�� ����������

% % % ��� ���� ������� �������� �������������� ���� isInTrace � SerNumber, 
% % % ��� ����: 
% % %    -1 ��������, ��� ������� ��������������� (������� ��� ����);
% % %     0 ��������, ��� ������ �� ������;
% % %     n > 0 ��������, ��� ������� ������ � ������ � ������� n, 
% % % � ����� ��� ������� ���������� ������ ������ ���� � ����� �������
% % % (� ������� ������ ��������� �������. ������ ��� ��������)
for i = 1:Nstep 
    Results(i).isInTrace = -1*ones(size(Results(i).Centers)); 
    %Results(i).SerNumber = zeros(size(Results(i).Centers)); 
end

% % �������� ������ ��� �������
% for i = 1:Nstep 
%     %Results(i).SerNumber = [Nplots+1 : Nplots+Results(i).Num]; 
%     Nplots = Nplots+Results(i).Num;
% end

%  ===== ��������������� ����� ���������������� ������� ===================
for i = 1:Nstep-2
    for k = 1:Results(i).Num % ������� ���� ������� (����������)
        %{
    ���� ����������� ���������� ������� ���������� � �����������
    � ��������� ���������?
    ������� �������� - ���� � �������� ����� ������� �������
     +/- deltas[1] �� ����� ���� ��������� �������, ��� ��������� ���, ��
     +/- deltas[2] �� ����� ���� ������� ����� ����
        %}
        Begn = Results(i).Begs(k); % ������ ������� �������
        Finh = Results(i).Fins(k); % ����� ������� �������
        Widh = Results(i).Widths(k); % ������������ ������� �������
        
        % ������� ��������� �� ���� ��� �����
        for k1 = 1:Results(i+1).Num
            Begn1 = Results(i+1).Begs(k1); % ������ ������� �������
            Finh1 = Results(i+1).Fins(k1); % ����� ������� �������
            Widh1 = Finh1-Begn1; % ������������ ������� �������
            % ������� ������������� � ������ ������� (deltas(1))?
            % ����.: ������, ������ ������, ������ ��������� �� ������� ���
            % � �������� ���
            if (Begn*(1-deltas(1)) > Finh1) || (Finh*(1+deltas(1)) < Begn1)
                % ������� ������� � ������� k1 �� �������� - ��������� � ���������
                continue
            else
                % ���� ��� ���������� - ������� ��������� � ������� �����
                Results(i).isInTrace(k) = 0;  break
            end
        end
        
        if (Results(i).isInTrace(k) == 0), continue, end % � ���������� k ��� i
        
        % ���� �� ������� ������ ����� ������� - ������� �� ��� ���� �����
        for k1 = 1:Results(i+2).Num
            Begn1 = Results(i+2).Begs(k1); % ������ ������� �������
            Finh1 = Results(i+2).Fins(k1); % ����� ������� �������
            Widh1 = Finh1-Begn1; % ������������ ������� �������
            % ������� ������������� � ������ ������� (deltas(2))?
            if (Begn*(1-deltas(2)) > Finh1) || (Finh*(1+deltas(2)) < Begn1)
                % ������� ������� � ������� k1 �� �������� - ��������� � ���������
                continue
            else
                % ���� ��� ���������� - ������� ��������� � ������� �����
                Results(i).isInTrace(k) = 0;  break
            end
        end
        
        % ���� � ����� ������� ���� �� ������� �������, �������
        % ������������� � �������, �� Results(i).isInTrace(k) �������
        % ������ -1 � ������� ������������� �� ����������� ������������
        
        % ���� �� �������� ����� ������������� ���������, � ������� ���-��
        % "�����" ���� �� ����������
    end
end

% ---- ������ ��� ������� - ������� ��������� ----------------------------
figure(999);
for k = 1:Nstep-2
    for tmp = 1:Results(k).Num
        if Results(k).isInTrace(tmp) >= 0
            center = round(Results(k).Centers(tmp));
            wdth   = Results(k).Widths(tmp);
            % ������� "��������" ������� ������� ������
            plot([k k], [Results(k).Begs(tmp) Results(k).Fins(tmp)],'r','LineWidth',2);
        end
    end
end
% ------------------------------------------------------------------------

% return

% ===== ���������� ������ ====================================
% ������� ������: ���� ��������� �������, � �������
% Results(i).isInTrace(k) = 0 � �������� �� �� �������� ������ ���, ��� ��
% ���� "�������" ������� ������ (��������� ������) ���� ����� �� ����� Ncptrd.
% ��� �������� ������ ��, ��� ������� �����������, �������� � ���-�� ����.
% ������� ��, ����� ������� ����������� ������ � ������ ������� ������� (������ ������������!)
% ���� ��������� � �������, ������� ��� ������ � �����-���� ������, �� (�
% ������� ������) �������� ������ ������������ � ��� ����� ����������� ��
% ��������� ������ ������� �������� ��� -1 ���� �� ������� ���� ���������� �����
% (���� ������������)
% ���� ��� ������ - �� ����� ������ �������� ������������� � ��������� ����� � ��
% ���������� � ������ ���� �� ������� �� ������� � ����� -1. 
% ������������� ������ ���������� ��������� ���� ���� �� �������� ��
% �������, ������� ��� ������ � ����� ���������� ������, ��� ���� ���
% ����� ������� isInTrace++. 
% �� ����, ��� ����� ��������� ��������, ������������� ��������� ������ ��
% ��-����� (������ ������ ������ ����� ���), ��, ��������, ��� "�����."
% ����� ����� ����� ������ ��������� ��������� �������� � ������������
% ����������� (����, ���� "�����" ������ ����� ������ ��-���������, ���
% ����� ���� �������� - ������������� ����� �������)
% ���� �� ������ ������� Ncptrd ������� �� ��������� ������, �� ���
% ������� �� ��������� ������ �������� -1 (����, ��������, ������ ��
% ���������� ����� � ��������� �������, ������� �������������� �����������
% ��������, �� ���� ������ ����� �� �����)

% ����������� ��� ��������� ������
traces = struct('Res_i',0,'Num_in_Res',0,'Center',0,'Clutter',0);
traces = repmat(traces, Nstep, MaxTracks); % ������� - ���� ������

Ntraces = 0; % ���-�� ������������ "������"

flag_too_many_traces = 0; % ��� ���������� ������ ��� ���������� ��������

% ������� �������, � ������� ����� ������ ��������� ������
for i = 1:Nstep-2
    for k = 1:Results(i).Num
        
        if (Results(i).isInTrace(k) == 0)
            
            % "��������" ��������� ������
            temp_trace  = struct('Res_i',0,'Num_in_Res',0,'Center',0,'Clutter',0);
            temp_trace  = repmat(temp_trace, Nstep, 1);
            temp_length = 1;
            isClutter   = 0;

            % ������ ��� ������� ��������� ������ ������� � ������
            %temp_lengs = zeros(Nstep, 1); 
            
            % ������ ������� - �������
            temp_trace(1).Res_i      = i;
            temp_trace(1).Num_in_Res = k;
            temp_trace(1).Center     = Results(i).Centers(k);
            %temp_lengs(1)            = Results(i).Widths(k);
            
            i2 = i; k2 = k; 
            
            % ������ �������� ���������� ��������� ������
            % ������� ���� ������������� � �� �������� � ������ ������ �������
            while(1)
                % ������� ��������� � ��������� ���������
                if ((temp_length > 1) && (isClutter == 1))
                    isClutter_prev = isClutter;
                else
                    isClutter_prev = 0;
                end
                i_prev = i2; k_prev = k2;
                [tmp,k2,isClutter] = GetNextClosest(Results(i_prev),k_prev,Results(i2+1:i2+2));
                i2 = i2+tmp;
                if (tmp > 0) % ���-������ �������?
                    %fprintf(1,'Temp trace #i %g, k = %g, k2 = %g, Results(i).isInTrace(k) = %g\n', i, k, k2, Results(i).isInTrace(k));
                    if (Results(i2).isInTrace(k2) > 0)
                        % ������� � �������, ������� ��� ����-�� ������ - �������� ��� ������� ��
                        % ��������� ������ -1 � ������ �� ����� ��������� �����
                        % �����, ������, ������������!
                        Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                        break
                        
                    elseif (temp_length >= Ncptrd)
                        % ������ ������������ - ��������� �� ��������� �� ���������� ���-�� �����
                        if (Ntraces >= MaxTracks)
                            warning('Too many traces have been detected!'); %#ok<WNTAG>
                            % ������� ����� ����� - ����� ��������
                            flag_too_many_traces = 1; 
                            break
                        else
                            % ����� - ������ ����� ������
                            Ntraces = Ntraces+1; % ����� ������
                            traces(:,Ntraces) = temp_trace;
                            Results = SetWholeTraceTo(Results,temp_trace,temp_length,Ntraces);
                            
                            % � ���������� ������ �� ���� ��������
                            while(1)
                                % ������� ��������� � ��������� ���������
                                isClutter_prev = isClutter; % �������� �������� �� ���� ���
                                if ((isClutter == 1))
                                    % ��������� ������� ������ � ��������� "������" - ����� � ������������
                                    i_prev = i2; % traces(tmp-1,Ntraces).Res_i; % ������ i2-1 ������ - ����� ������� ��������
                                    k_prev = k2; %traces(tmp,Ntraces).Num_in_Res; % �� �� �����
                                    %i2 = i2+1; % ������������� ����� ���������� ������
                                else
                                    i_prev = i2; k_prev = k2;
                                end
                                [tmp,k2,isClutter] = GetNextClosest(Results(i_prev),k_prev,Results(i2+1:i2+2));
                                i2 = i2+tmp; % �� ������� ����� ��������� (1 ��� 2, 2 ��� 3, ���� ���������� ��� �������)
                                if (tmp > 0)
                                    % ���� ������ ��������� (�������) �����
                                    tmp = find(~[traces(:,Ntraces).Res_i],1); % �������� ������
                                    traces(tmp,Ntraces).Res_i       = i2;
                                    traces(tmp,Ntraces).Num_in_Res  = k2;
                                    traces(tmp,Ntraces).Center      = Results(i2).Centers(k2); % ����� ����� ����������, ���� isClutter (��. �����)
                                    traces(tmp,Ntraces).Clutter     = isClutter; % ������� "������"
                                    Results(i2).isInTrace(k2)       = Ntraces;
                                    % ����������� ���������� �������, ���� ��� ������ � ��������� "������", ����� ������ ����� ���� � ������� 
                                    if (isClutter),
                                        x = [traces(tmp-Nlin : tmp-1, Ntraces).Res_i];
                                        Y = [traces(tmp-Nlin : tmp-1, Ntraces).Center];
                                        xi = traces(tmp, Ntraces).Res_i;
                                        % ������ ������������ ������ ����������� �������� ������������� �� Nlin ���������� ��������
                                        traces(tmp,Ntraces).Center = interp1(x, Y, xi,'linear','extrap');
                                    end
                                    
                                    % ---- ������ ��� ������� --------------------------------------
                                    fprintf(1,'================================\n');
                                    fprintf(1,'x = %g\n', traces(tmp,Ntraces).Res_i);
                                    fprintf(1,'y = %g\n', traces(tmp,Ntraces).Num_in_Res);
                                    fprintf(1,'Center = %g\n', traces(tmp,Ntraces).Center);
                                    fprintf(1,'isClutter %g\n', traces(tmp,Ntraces).Clutter);
                                    fprintf(1,'Trace No.%g\n', Ntraces);
                                    % ------------------------------------------------------------------------
                                    
                                    % ��������� ������ ���� ���������� ������
                                    temp_length                        = temp_length+1; % ��������� ����� ������
                                    temp_trace(temp_length).Res_i      = i2;
                                    temp_trace(temp_length).Num_in_Res = k2;
                                    temp_trace(temp_length).Center     = traces(tmp,Ntraces).Center; % � �����!!!
                                    temp_trace(temp_length).Clutter    = isClutter;
                                    %temp_lengs(temp_length)            = Results(i2).Widths(k2);

                                    if (i2 >= length(Results)-3), % ��� �� �� �������� �� ������� ������������� 
                                        ShowMeTheTrace(Results,traces,Ntraces);
                                        break; 
                                    end
                                else
                                    % ����� ���������� - ������ ������ � ��������� � ������ ����� ��������� �����
                                    %Results = SetWholeTraceTo(Results,traces(:,Ntraces),tmp,Ntraces);
                                    ShowMeTheTrace(Results,traces,Ntraces);
                                    %pause;
                                    break % ������ � ������� Ntraces ���������
                                end
                            end
                            
                            % ������ �� ����� ��������� ������
                            break
                        end
                        
                    else
                        % ���������� ����������� ��������� ������
                        temp_length                        = temp_length+1;
                        temp_trace(temp_length).Res_i      = i2;
                        temp_trace(temp_length).Num_in_Res = k2;
                        temp_trace(temp_length).Clutter    = isClutter;
                        %temp_lengs(temp_length)              = Results(i2).Widths(k2);
                        temp_trace(temp_length).Center     = Results(i2).Centers(k2); % � ��� ����!
                        
                        if (i2 >= length(Results)-2),
                            % ��������� ������ �� ������� - �������������
                            % ���������
                            Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                            break;
                        end
                        
                    end
                else
                    % ��� �������� ����������� - �������� ��� ������� ��
                    % ��������� ������ -1 � ������ �� ����� ���������
                    % �����
                    Results = SetWholeTraceTo(Results,temp_trace,temp_length,-1);
                    %flag_to_exit = 1; 
                    break
                end
                
            end % of while
            
            if (flag_too_many_traces), break, end
        end
        % ����� ����� �������� ��������� ������
        
        
    end
    if (flag_too_many_traces), break, end
end

% ---- ������ ��� ������� - ������� ��������� ----------------------------
% figure(999);
% for i = 1:Ntraces,
%     ShowMeTheTrace(Results, traces, i);%
% end
% figure(999); xlabel('Classified: RAMERA type radar','FontSize', 16);
% ------------------------------------------------------------------------

% ==== ��������� �������� ��������� ������ ===============================
% ���������: ������� �������, ������������ ��������� �������, �������
% "�������", �������� � ������������ ����� �������������� �������������
% ����� (��������� �� ���������� �� ���), ��� �������, ��������� ���
% ������� ����� � ����� ����� (���������� �� ������-"�����������"),
% ����� �����, etc.

% ��-������, ��������� ����� ������������ ������, ���� �� ������� �����,
% �� ���� ��� �������� ����� ������, ���� ������� ����� ��� ���
% �������������. ����� ������ - ���������������.
if (Ntraces > MaxTraces)
    fprintf(1,'===== Final Result ====================================\n');
    fprintf(1,'Too many traces: unclassified.\n');
    return
end

% ���� ����� �������� ��������� ������
TracksInfo = struct('Fmean',0,'dFmax',0,'Length',0,'Curvature',0,'MSE',0,'isBreak',0);
TracksInfo = repmat(TracksInfo, Ntraces);
Lengths    = zeros(Ntraces, 1); % ����� �����

for n_trk = 1:Ntraces
    the_last = find(~[traces(:,n_trk).Res_i],1)-1; % ����� ������� ������
    Lengths(n_trk) = the_last;
    
    x        = [traces(1:the_last,n_trk).Res_i];
    tmp      = [traces(1:the_last,n_trk).Num_in_Res];
    wtdhs    = zeros(size(x));
    
    for i = 1:length(x)
        wtdhs(i) = Results(x(i)).Types(tmp(i)); % ���� ������� �� ����������������
    end
    
    % �������� �� ������� - ������������ ���� � ���� ��-�����
    if (mean(wtdhs)/length(wtdhs) > 0.5),
        fprintf(1,'===== Final Result ====================================\n');
        fprintf(1,'Most probable radar: %s\n', 'Stalker 34G');
        fprintf(1,'With probability: %g%%\n', mean(wtdhs)/length(wtdhs));
        return % ������� - ������ ����������� �������, ������ �� � ���� ���� ���
    end
    
    % �������� ���� �������
    y = [traces(1:the_last,n_trk).Center]; % ������ �������
    
    % �������� �������� ���������
    TracksInfo(n_trk).Fmean  = mean(y);
    TracksInfo(n_trk).dFmax  = max(y)-min(y);
    TracksInfo(n_trk).Length = x(end)-x(1);
    p = polyfit(x,y,2); TracksInfo(n_trk).Curvature = p(1)/p(2);
    p = polyfit(x,y,10); Y = polyval(p,x); err = abs(Y-y);
    TracksInfo(n_trk).MSE = std(err);
    
    % �������� �� ������� "�������" ������
    trh_brk = mean(err)+3.5*std(err); % ����� ����������� �������
    
    if ~isempty(find(err > trh_brk,1))
        % ���� � ������� - ����� ����������������� �������
        TracksInfo(n_trk).isBreak = 1; % ��������� ����� �����
        [~,i_max] = max(err); % �������� ��������� ������ �����
    end
    
end

% ---- ������ ��� ������� - ������� ��������� ������������� --------------
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

% ---- ���������� ������������� � ����� ����������� ----------------------
if (Ntraces == 1)
    if (TracksInfo(1).dFmax < 10)
        % ���-�� ������������ (�������� "���������� 6F", �� ���)
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
    % ���� ������ ������, �� ������� ������ ��������� ��� �������������
    if (IntersectEnough(traces(:,1:Ntraces),Lengths))
        % ���� ��� ���������� ������� � ������ ������������� - �� ���,
        % ������ ����� "������-CD"
        Type_of_radar = 'Multa CD';
        tmp = [TracksInfo(:).Length];
        Prob_of_recogn = min(tmp(1:2))/max(tmp(1:2));
    else
        if ((Ntraces == 2) || (Ntraces == 3))
            % ���� ����� �� ������������� � ������ ��� ��� ���, �� ���,
            % ������ ����� "������", �� ������� �� ���� ��� ������ �����
            % ���� �� �����
            Type_of_radar = 'Iskra';
            tmp = [TracksInfo(:).Length];
            Prob_of_recogn = (traces(1,2).Res_i-traces(Lengths(1),1).Res_i)/Nstep;
        else
            % ���-�� �� ��������
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



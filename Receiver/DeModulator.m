function [data_DeCoded_Rx,data_Rx,t2_rx ] = DeModulator(Symbol_stream_Rx,Coding_config,OFDM_config,Simulation_config,Audio_config )
%DEMODULATOR Summary of this function goes here
%   Detailed explanation goes here

%% Inputs
Coding_flag=Coding_config.Coding_flag;
K_coding=Coding_config.K_coding;
N_coding=Coding_config.N_coding;
Interleave_flag=Coding_config.Interleave_flag;

M=OFDM_config.M;
N_data=OFDM_config.N_data;
P_data=OFDM_config.P_data;
%N_preamble_CE=OFDM_config.N_preamble_CE;
%N_preamble_synch=OFDM_config.N_preamble_synch;



Voice_flag=Audio_config.Voice_flag;


N_symbols=Simulation_config.N_symbols;
Identical_symbols=Simulation_config.Identical_symbols;

if ~isempty(Simulation_config.Voice_only)
    Voice_only=Simulation_config.Voice_only;
else
    Voice_only=0;
end

%%

hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'NormalizationMethod','Average Power','AveragePower',P_data,'BitOutput',1);

Symbol_stream_Rx=Symbol_stream_Rx*sqrt(P_data/var(Symbol_stream_Rx)); % normlization to P_data. needed in 2 PC's link
data_Rx=step(hDemod,Symbol_stream_Rx);

%%
if Voice_flag
    q_bits=Audio_config.q_bits;
    Fs_Audio=Audio_config.Fs_Audio;
    
    q=quantizer([q_bits q_bits-1]);
    Speaker=dsp.AudioPlayer('SampleRate',Fs_Audio); % must have the same Fs as Microphone so that the voice does not distort
    
    WAV_read=dsp.AudioFileReader('Filename','Lior_recovered.wav'...
        ,'OutputDataType','int16');
else
    N_samples_frame=[];
end


if Coding_flag
    
    if Interleave_flag
        data_Rx=matdeintrlv(data_Rx,length(data_Rx)/N_coding,N_coding); % matrix interleaving; accumulates the whole signal into a matrix row by row and then transfres it on column by column. the receive packet needs to have the same size of the transmit packet (bit-wise or symbol-wise) so that the matrix de-intrleaving operates correctly
    end
    
    %%% Decode
    dec = comm.RSDecoder('BitInput',true,'CodewordLength',N_coding,'MessageLength',K_coding);
    data_DeCoded_Rx = step(dec,data_Rx); % RS decode
    
    if Voice_flag
        N_bits=ceil(N_symbols/(K_coding*N_data*q_bits))*(N_data*K_coding*q_bits)*log2(M);
        
        data_DeCoded_Rx2=reshape(data_DeCoded_Rx,8,[]);
        data_DeCoded_Rx1=data_DeCoded_Rx2';
        Signal_voice_compressed=bi2de(data_DeCoded_Rx1);
        
        fwrite(fopen('Lior_r_w.bin','w'),Signal_voice_compressed,'uint8')
        
        [a,b]=system('cp_amrwb_decoder Lior_r_w.bin Lior_recovered.wav');
        
        t2_rx=clock; % communication terminates here
        
        B=0;
        while ~isDone(WAV_read)
            %     b=step(WAV_read);
            %     B=[B;b];
            B=step(WAV_read);
            step(Speaker,B);
        end
        release(Speaker)
        
        if Voice_only
            stop
        end
        
    else % Not voice
        
        t2_rx=clock; % communication terminates here
    end
    
    
%     %%% Error calculations
%     hError=comm.ErrorRate;
%     Error_strct=step(hError,data_PreCoded_Tx,data_DeCoded_Rx); % we ommit from the comparison the chopped symbols (due to group delay)
    
else % Not Coding

    if Voice_flag

        
        data_Rx2=reshape(data_Rx,8,[]);
        data_Rx1=data_Rx2';
        Signal_voice_compressed=bi2de(data_Rx1);
        
        fwrite(fopen('Lior_r_w.bin','w'),Signal_voice_compressed,'uint8')
        
        [a,b]=system('cp_amrwb_decoder Lior_r_w.bin Lior_recovered.wav');
        
        t2_rx=clock; % communication terminates here
        
        B=0;
        while ~isDone(WAV_read)
            %     b=step(WAV_read);
            %     B=[B;b];
            B=step(WAV_read);
            step(Speaker,B);
        end
        release(Speaker)
        
        if Voice_only
            quit
        end

    else
        t2_rx=clock; % communication terminates here
        
        
    end
 
    data_DeCoded_Rx=[];
%     
%     hError=comm.ErrorRate;
%     Error_strct=step(hError,data_Tx,data_Rx); % we ommit from the comparison the chopped symbols: due to group delay an fgthe fact that matlab's filter function chops the end of the vector, and we connot do anything about that
    
end




end


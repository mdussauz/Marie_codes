Across_mice_waveform_correlation_script

%% Classic open loop sessions
SessionName = {'O3/O3_20211005_r0_processed.mat',...
'O8/O8_20220702_r0_processed.mat',...
'O9/O9_20220630_r0_processed.mat',...
'S1/S1_20230314_r0_processed.mat',...
'S3/S3_20230321_r0_processed.mat',...
'S6/S6_20230718_r0_processed.mat',... % !!!this session is a free lever one - to be changed when classic is fixed!!!
'S7/S7_20230707_r0_processed.mat',... 
'S11/S11_20230812_r0_processed.mat',...
'S12/S12_20230727_r0_processed.mat'};

[WFcorrelation, concWFcorrelation, CorrelationCriteria, concCorrelationCriteria] = GetWaveformCorrelation (SessionName, sessionHalfTime); 

nb_units_total(mouse) = length(concWFcorrelation);
%based on concatenated waveforn across the channels of the unit's main
%tetrode:
bad_units_conc = find(concWFcorrelation > concCorrelationCriteria);
nb_units_conc(mouse) = sum(concWFcorrelation > concCorrelationCriteria);

%based on single waveform on the unit's main channel:
bad_units_single = find(WFcorrelation > CorrelationCriteria);
nb_units_single(mouse) = find(WFcorrelation > CorrelationCriteria);


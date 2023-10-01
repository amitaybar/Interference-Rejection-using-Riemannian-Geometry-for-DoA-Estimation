
%% Plot Params
linewd                  = 0.8;
hcfontsize              = 20;
MarkerSize              = 9;
%%
if PlogFigs_3_4_10
    %%
    SirMeanMLAll = [SirMeanMLRMScen1SNR20(:) SirMeanMLEMScen1SNR20(:) SirMeanMLRMScen12SNR20(:) SirMeanMLEMScen12SNR20(:) ...
        SirMeanMLRMScen13SNR20(:) SirMeanMLEMScen13SNR20(:)];
    DirectivityMLAll = [DirectivityMLRSaveMScen1SNR20(:) DirectivityMLESaveMScen1SNR20(:) DirectivityMLRSaveMScen12SNR20(:) DirectivityMLESaveMScen12SNR20(:) ...
        DirectivityMLRSaveMScen13SNR20(:) DirectivityMLESaveMScen13SNR20(:)];
    
    
    LabelsML1 = {'Riem.','Euc.', 'Riem.' , 'Euc.', 'Riem.' , 'Euc.'};
    LabelsML2 = {'-6dB','-6dB', '-10dB','-10dB','-20dB','-20dB' };
    LabelsAllML = {LabelsML1,LabelsML2};
    save(['DirectivityMLAll_Mscen1'],'DirectivityMLAll')
    save(['LabelsAllML_Mscen1'],'LabelsAllML')
    MyBoxSubPlotFunc(DirectivityMLAll, LabelsAllML, 'Directivity - DS' ,3);
    saveas(gcf,fullfile('results', 'Fig3Bottom_DirectivityDS.jpg'))    
    save(['SirMeanMLAll_Mscen1'],'SirMeanMLAll')
    save(['LabelsAllML_Mscen1'],'LabelsAllML')
    MyBoxSubPlotFunc(SirMeanMLAll, LabelsAllML, 'SIR - DS' ,3);
    saveas(gcf,fullfile('results', 'Fig3Top_SirDS.jpg'))    
    %
    %-------------------------------------------------------------------------------------
    DirectivitySSAll = [DirectivityMusicRSaveMScen1SNR20(:) DirectivityMusicESaveMScen1SNR20(:) DirectivityIntersectSaveMScen1SNR20(:)....
        DirectivityMusicRSaveMScen12SNR20(:) DirectivityMusicESaveMScen12SNR20(:)    DirectivityIntersectSaveMScen12SNR20(:) ...
        DirectivityMusicRSaveMScen13SNR20(:) DirectivityMusicESaveMScen13SNR20(:)    DirectivityIntersectSaveMScen13SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Intersection', 'Riem.' , 'Euc.', 'Intersection', 'Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'-6dB','-6dB', '-6dB', '-10dB','-10dB', '-10dB','-20dB','-20dB' ,'-20dB'};
    LabelsAllSS = {LabelsML1,LabelsML2};  
    
    
    %-------------------------------------------------------------------------------------
    DirectivityOracleSSAll = [DirectivityMusicRSaveOracleMScen1SNR20(:) DirectivityMusicESaveOracleMScen1SNR20(:) DirectivityIntersectSaveOracleMScen1SNR20(:)....
        DirectivityMusicRSaveOracleMScen12SNR20(:) DirectivityMusicESaveOracleMScen12SNR20(:)    DirectivityIntersectSaveOracleMScen12SNR20(:) ...
        DirectivityMusicRSaveOracleMScen13SNR20(:) DirectivityMusicESaveOracleMScen13SNR20(:)    DirectivityIntersectSaveOracleMScen13SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Intersection', 'Riem.' , 'Euc.', 'Intersection', 'Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'Oracle','Oracle','Oracle', 'Oracle' , 'Oracle', 'Oracle', 'Oracle' , 'Oracle', 'Oracle'};
    LabelsML3 = {'-6dB','-6dB', '-6dB', '-10dB','-10dB', '-10dB','-20dB','-20dB' ,'-20dB'};
    LabelsAllSS = {LabelsML1,LabelsML2,LabelsML3};
    % 
    %---------------------------------------------------------------------------------------
    SirMeanSSAll = [SirMeanMusicRMScen1SNR20(:) SirMeanMusicEMScen1SNR20(:) SirMeanIntersectMScen1SNR20(:)....
        SirMeanMusicRMScen12SNR20(:) SirMeanMusicEMScen12SNR20(:)    SirMeanIntersectMScen12SNR20(:) ...
        SirMeanMusicRMScen13SNR20(:) SirMeanMusicEMScen13SNR20(:)    SirMeanIntersectMScen13SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Intersection', 'Riem.' , 'Euc.', 'Intersection', 'Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'-6[dB]','-6[dB]', '-6[dB]', '-10[dB]','-10[dB]', '-10[dB]','-20[dB]','-20[dB]' ,'-20[dB]'};
    LabelsAllSS = {LabelsML1,LabelsML2};
    save(['SirMeanSSAll_Mscen1_12_13'],'SirMeanSSAll')
    save(['LabelsAllSS_Mscen1_12_13'],'LabelsAllSS')
    MyBoxSubPlotFunc(SirMeanSSAll, LabelsAllSS, 'SIR - SbSp ' ,3);
    saveas(gcf,fullfile('results', 'Fig4Top_SirSbSp.jpg'))    
    %-----------------------------------------------------------------------------------------
    SirMeanSSORacleAll = [SirMeanMusicROracleMScen1SNR20(:) SirMeanMusicEOracleMScen1SNR20(:) SirMeanIntersectOracleMScen1SNR20(:)....
        SirMeanMusicROracleMScen12SNR20(:) SirMeanMusicEOracleMScen12SNR20(:)    SirMeanIntersectOracleMScen12SNR20(:) ...
        SirMeanMusicROracleMScen13SNR20(:) SirMeanMusicEOracleMScen13SNR20(:)    SirMeanIntersectOracleMScen13SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Intersection', 'Riem.' , 'Euc.', 'Intersection', 'Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'Oracle','Oracle','Oracle', 'Oracle' , 'Oracle', 'Oracle', 'Oracle' , 'Oracle', 'Oracle'};
    LabelsML3 = {'-6dB','-6dB', '-6dB', '-10dB','-10dB', '-10dB','-20dB','-20dB' ,'-20dB'};
    LabelsAllSS = {LabelsML1,LabelsML2,LabelsML3};
    save(['SirMeanSSORacleAll_Mscen1_12_13'],'SirMeanSSORacleAll')
    save(['LabelsAllSS_Oracle_Mscen1_12_13'],'LabelsAllSS')
    MyBoxSubPlotFunc(SirMeanSSORacleAll, LabelsAllSS, 'SIR - SbSp Oracle' ,3);
    saveas(gcf,fullfile('results', 'Fig4Bottom_SirSbSpOracle.jpg'))    
    %------------------------------------------------------------------------------------------------------------------
    DirectivityMvdrAll = [DirectivityMvdrRSaveMScen1SNR20(:) DirectivityMvdrESaveMScen1SNR20(:) DirectivityMvdrRSaveMScen12SNR20(:) DirectivityMvdrESaveMScen12SNR20(:) ];
    LabelsMvdr1 = {'Riem.','Euc.', 'Riem.' , 'Euc.', 'Riem.' , 'Euc.'};
    LabelsMvdr11 = {'Riem.','Euc.'};%, 'Riem.' , 'Euc.'};%, 'Riem.' , 'Euc.'};
    LabelsMvdr2 = {'-6dB','-6dB', '-10dB','-10dB','-20dB','-20dB' };
    LabelsAllMvdr  = {LabelsMvdr1,LabelsMvdr2};
    % 
    
    SirMvdrAll = [SirMeanMvdrRMScen1SNR20(:) SirMeanMvdrEMScen1SNR20(:) SirMeanMvdrRMScen12SNR20(:) SirMeanMvdrEMScen12SNR20(:) SirMeanMvdrRMScen13SNR20(:) SirMeanMvdrEMScen13SNR20(:)];
    LabelsMvdr1 = {'Riem.','Euc.', 'Riem.' , 'Euc.', 'Riem.' , 'Euc.'};
    LabelsMvdr11 = {'Riem.','Euc.'};%, 'Riem.' , 'Euc.'};%, 'Riem.' , 'Euc.'};
    LabelsMvdr2 = {'-6dB','-6dB', '-10dB','-10dB','-20dB','-20dB' };
    LabelsAllMvdr  = {LabelsMvdr1,LabelsMvdr2};
    save(['SirMvdrAll_Mscen1'],'SirMvdrAll')
    save(['LabelsAllMvdr_Mscen1'],'LabelsAllMvdr') 
    MyBoxSubPlotFunc(SirMvdrAll, LabelsAllMvdr, 'SIR - MVDR' ,3,[2 2 2]);
    saveas(gcf,fullfile('results', 'Fig9_MVDR.jpg'))    
end
if PlotFigs_8
    %-------------------------------------------------------------------------------------------------------------
    SirMeanAll = [SirMeanMLRMScen27SNR20(:) SirMeanMLEMScen27SNR20(:) SirMeanMusicRMScen27SNR20(:) SirMeanMusicEMScen27SNR20(:) SirMeanIntersectMScen27SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'DS','DS','SbSp' , 'SbSp', 'SbSp'};
    LabelsAllML = {LabelsML1,LabelsML2};
    DirectivityAll = [DirectivityMLRSaveMScen27SNR20(:) DirectivityMLESaveMScen27SNR20(:) DirectivityMusicRSaveMScen27SNR20(:) DirectivityMusicESaveMScen27SNR20(:) DirectivityIntersectSaveMScen27SNR20(:)];
    
    save(['SirMeanAll_Mscen27'],'SirMeanAll')
    save(['LabelsAllML_Mscen27'],'LabelsAllML')
    MyBoxSubPlotFunc(SirMeanAll, LabelsAllML, 'SIR - Multiple interferences ' ,2,[2,3]);
    saveas(gcf,fullfile('results', 'Fig7Right_Sir_14Interferences.jpg'))
    % 
end
%---------------------------------------------------------------------------------------------------------------
if exist('SirMeanMLRMScen29SNR20','var')
    SirMeanAll = [SirMeanMLRMScen29SNR20(:) SirMeanMLEMScen29SNR20(:) SirMeanMusicRMScen29SNR20(:) SirMeanMusicEMScen29SNR20(:) SirMeanIntersectMScen29SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'DS','DS','SbSp' , 'SbSp', 'SbSp'};
    LabelsAllML = {LabelsML1,LabelsML2};
    DirectivityAll = [DirectivityMLRSaveMScen29SNR20(:) DirectivityMLESaveMScen29SNR20(:) DirectivityMusicRSaveMScen29SNR20(:) DirectivityMusicESaveMScen29SNR20(:) DirectivityIntersectSaveMScen29SNR20(:)];
    
    MyBoxSubPlotFunc(SirMeanAll, LabelsAllML, 'SIR - Multiple interferences ' ,2,[2,3]);
    
    MyBoxSubPlotFunc(DirectivityAll, LabelsAllML, 'Directivity - Multiple interferences' ,2,[2 3]);
end
%---------------------------------------------------------------------------------------------------------------
if PlotFigs_9
    SirMeanAll = [SirMeanMLRMScen20SNR20(:) SirMeanMLEMScen20SNR20(:) SirMeanMusicRMScen20SNR20(:) SirMeanMusicEMScen20SNR20(:) SirMeanIntersectMScen20SNR20(:) ];
    DirectivityAll = [DirectivityMLRSaveMScen20SNR20(:) DirectivityMLESaveMScen20SNR20(:) DirectivityMusicRSaveMScen20SNR20(:) DirectivityMusicESaveMScen20SNR20(:) DirectivityIntersectSaveMScen20SNR20(:)];
    LabelsML1 = {'Riem.','Euc.','Riem.' , 'Euc.', 'Intersection'};
    LabelsML2 = {'DS','DS','SbSp' , 'SbSp', 'SbSp'};
    LabelsAllML = {LabelsML1,LabelsML2};
    save(['SirMeanAll_Mscen20'],'SirMeanAll')
    save(['LabelsAllML_Mscen20'],'LabelsAllML')
    MyBoxSubPlotFunc(SirMeanAll, LabelsAllML, 'SIR - Two targets' ,2,[2 3]);
    saveas(gcf,fullfile('results', 'Fig8Right_Sir_2desired.jpg'))   
    %------------------------------------------------------------------------------------------------------------------
    SirMeanAll = [SirMeanMLRMScen21SNR20(:) SirMeanMLEMScen21SNR20(:) SirMeanMusicRMScen21SNR20(:) SirMeanMusicEMScen21SNR20(:) SirMeanIntersectMScen21SNR20(:) ];
    DirectivityAll = [DirectivityMLRSaveMScen21SNR20(:) DirectivityMLESaveMScen21SNR20(:) DirectivityMusicRSaveMScen21SNR20(:) DirectivityMusicESaveMScen21SNR20(:) DirectivityIntersectSaveMScen21SNR20(:)];
    save(['SirMeanAll_Mscen21'],'SirMeanAll')
    save(['LabelsAllML_Mscen21'],'LabelsAllML')
    MyBoxSubPlotFunc(SirMeanAll, LabelsAllML, 'SIR - Two targets' ,2,[2 3]);
    saveas(gcf,fullfile('results', 'Fig8Left_Sir_2desired.jpg'))
    % 
end

%---------------------------------------------------------------------------------------------------------------
%
function [daily] = daily(Data)
%  DAILY FUNCTION ANALYZING BEHAVIORAL RHYTHMS OF PERFORMANCE - DAILY
% SOME GLOBAL VARIABLES
minRT           = 150;
maxRT           = 800;
%dailyPlots      = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\DAILY';
dailyPlots      = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\DAILY';
% ASSIGNED TO FREDDIES PLOTS
%dailyPlots      = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\FreddieDAILY';
%dailyPlots      = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\FreddieDAILY';
% READ FILES TO SEE WHAT TO SKIP
% SINGLE SESSION RHYTHMIC MEASURES
sessions = unique(Data.date);
CTDmin      = 500;
CTDmax      = 1400;
CTDstep     = 10;
Binwidth    = 50;
Bins        = ((CTDmin+(CTDmin+Binwidth))/2):CTDstep:((CTDmax+(CTDmax-Binwidth))/2);
NumBins     = length(Bins);
for i = 1:length(sessions)
    DailyTbl                = Data(Data.date == sessions(i),:);                 % USING FOR INDIVIDUAL SESSION RHYTHMIC METRICS & PLOTS
    EyeBrkIdx               = ismember(DailyTbl.break_eye,'true');              % REMOVING EYE BREAKS FROM THE START
    PreCueLevBrkIdx         = ismember(DailyTbl.break_eye,'false') & ismember(DailyTbl.break_lever,'true') & (DailyTbl.t_cue_on == 0);  % LEVER BREAK TRIALS BEFORE A CUE ONSET; INDEPENDENT OF EyeBrkIdx 
    RmvIdx                  = logical(EyeBrkIdx + PreCueLevBrkIdx);
    UseableTbl              = DailyTbl;                                         % MAKING SEPARATE TABLE INSTEAD OF MODIFYING DailyTbl. I WANT TO KEEP DAILY TABLE INTACT, SO THAT I CAN TRACK THOSE REMOVE TRIALS LATER ON IF NEEDED
    UseableTbl(RmvIdx,:)    = [];
    NumTrls                 = height(UseableTbl);
    % TABLE SUBDIVISIONS (VALID, INVALID S.O., INVALID D.O.)
    % VALID
    ValIdx                  = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 1:8;    %ISOLATE VALID CUE TRIALS
    ValIdx                  = logical(sum(ValIdx,2));
    ValTbl                  = UseableTbl(ValIdx,:);     % MAKE VALID CUED TABLE
    ValCorrIdx              = ((ValTbl.t_contrast_change + minRT) <= ValTbl.t_release) & (ValTbl.t_release <= (ValTbl.t_contrast_change + maxRT));
    ValCorrTbl              = ValTbl(ValCorrIdx,:);     % ISOLATE VALID CUE TRIALS THAT ARE CORRECT
    ValIncorrTbl            = ValTbl(ValCorrIdx==0,:);  % ISOLATE VALID CUE TRIALS THAT ARE INCORRECT
    ValIncorrTbl            = ValIncorrTbl(ismember(ValIncorrTbl.break_lever,'false'),:);   % REMOVES LEVER BREAK TRIALS (CAN'T BE CALCULATED)
    % INVALID SAME OBJECT (S.O.)
    InvalIdxSO              = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 9;  %ISOLATE INVALID S.O. TRIALS
    InvalTblSO              = UseableTbl(InvalIdxSO,:); % MAKE INVALID S.O. TABLE
    InvalCorrIdxSO          = ((InvalTblSO.t_contrast_change + minRT) <= InvalTblSO.t_release) & (InvalTblSO.t_release <= (InvalTblSO.t_contrast_change + maxRT)); 
    InvalCorrTblSO          = InvalTblSO(InvalCorrIdxSO,:);
    InvalIncorrTblSO        = InvalTblSO(InvalCorrIdxSO == 0,:);
    InvalIncorrTblSO        = InvalIncorrTblSO(ismember(InvalIncorrTblSO.break_lever,'false'),:);   % REMOVES LEVER BREAK TRIALS (CAN'T BE CALCULATED)
    % INVALID DIFFERENT OBJECT (D.O.)
    InvalIdxDO              = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 10; %ISOLATE INVALID D.O. TRIALS
    InvalTblDO              = UseableTbl(InvalIdxDO,:); % MAKE INVALID D.O. TABLE
    InvalCorrIdxDO          = (InvalTblDO.t_contrast_change + minRT <= InvalTblDO.t_release) & (InvalTblDO.t_release <= (InvalTblDO.t_contrast_change + maxRT));
    InvalCorrTblDO          = InvalTblDO(InvalCorrIdxDO,:);
    InvalIncorrTblDO        = InvalTblDO(InvalCorrIdxDO == 0,:);
    InvalIncorrTblDO        = InvalIncorrTblDO(ismember(InvalIncorrTblDO.break_lever,'false'),:);   % REMOVES LEVER BREAK TRIALS (CAN'T BE CALCULATED)
    % GENERATE CUE TARGET INTERVALS FOR THE RESPECTIVE TABLES
    ValCorrCueTar           = ValCorrTbl.t_contrast_change - ValCorrTbl.t_cue_on;
    ValIncorrCueTar         = ValIncorrTbl.t_contrast_change - ValIncorrTbl.t_cue_on;
    InvCorrSO_CueTar        = InvalCorrTblSO.t_contrast_change - InvalCorrTblSO.t_cue_on;
    InvIncorrSO_CueTar      = InvalIncorrTblSO.t_contrast_change - InvalIncorrTblSO.t_cue_on;
    InvCorrDO_CueTar        = InvalCorrTblDO.t_contrast_change - InvalCorrTblDO.t_cue_on;
    InvIncorrDO_CueTar      = InvalIncorrTblDO.t_contrast_change - InvalIncorrTblDO.t_cue_on;
    % LOOP FOR HIT RATES BASED ON CUE-TARG INTERVAL RANGES
    strtTme                 = CTDmin;  % FIRST CUE-TARG INT BEGINNING TIME %CHANGES DEPENDING ON WHAT'S BEING LOOKED AT
    endTme                  = CTDmin+Binwidth;  % FIRST CUE-TARG INT END TIME %CHANGES DEPENDING ON WHAT'S BEING LOOKED AT
    slide                   = CTDstep;   % VALUE IN ms TO MOVE THE SLIDING WINDOW FORWARD
    maxTime                 = CTDmax; % MAXIMUM CUE-TARG INT
    %dur                     = (strtTme:endTme);
    ValHR_Rthm              = zeros(1,length(endTme:slide:maxTime));
    InvSO_HR_Rthm           = zeros(1,length(endTme:slide:maxTime));
    InvDO_HR_Rthm           = zeros(1,length(endTme:slide:maxTime));
    CtrTmeWindow            = zeros(1,length(endTme:slide:maxTime));
    iter                    = 1;
    while endTme <= maxTime
        ValCorr             = sum(ValCorrCueTar >= strtTme & ValCorrCueTar <= endTme);
        ValIncorr           = sum(ValIncorrCueTar >= strtTme & ValIncorrCueTar <=endTme);
        ValHR               = ValCorr/(ValCorr + ValIncorr);
        ValHR_Rthm(iter)    = ValHR;
        InvCorrSO           = sum(InvCorrSO_CueTar >= strtTme & InvCorrSO_CueTar <= endTme);
        InvIncorrSO         = sum(InvIncorrSO_CueTar >= strtTme & InvIncorrSO_CueTar <= endTme);
        InvHRSO             = InvCorrSO/(InvCorrSO + InvIncorrSO);
        InvSO_HR_Rthm(iter) = InvHRSO;
        InvCorrDO           = sum(InvCorrDO_CueTar >= strtTme & InvCorrDO_CueTar <= endTme);
        InvIncorrDO         = sum(InvIncorrDO_CueTar >= strtTme & InvIncorrDO_CueTar <= endTme);
        InvHRDO             = InvIncorrDO/(InvCorrDO + InvIncorrDO);
        InvDO_HR_Rthm(iter) = InvHRDO;
        
        CtrTmeWindow(iter)  = (strtTme+endTme)/2;
        strtTme             = strtTme + slide;
        endTme              = endTme + slide;
        iter                = iter+1;
    end
    
    % TRENDED
    ValHR_Rthm(isnan(ValHR_Rthm))       = 0;
    InvSO_HR_Rthm(isnan(InvSO_HR_Rthm)) = 0;
    InvDO_HR_Rthm(isnan(InvDO_HR_Rthm)) = 0;
    %DETRENDED
    ValHR_RthmDtnd          = detrend(ValHR_Rthm);
    InvSO_HR_RthmDtnd       = detrend(InvSO_HR_Rthm);
    InvDO_HR_RthmDtnd       = detrend(InvDO_HR_Rthm);
    % DATE STRING
    dateStr                 = datetime(sessions(i),'InputFormat','yyMMdd');
    dateStr                 = string(dateStr);
    % Subplots ('Tiled Layout')
    figure('Position', [1,1,2100,600])
    t = tiledlayout(1,3);
    
    %TILE 1
    nexttile
    line(1)                 = plot(Bins,InvDO_HR_Rthm,'LineWidth',1.1,'Color',"#fdc086",'DisplayName','Invalid DO Rhythm');
    hold on
    line(2)                 = plot(Bins,InvSO_HR_Rthm,'LineWidth',1.1,'Color',"#7fc97f",'DisplayName','Invalid SO Rhythm');
    hold on
    line(3)                 = plot(Bins,ValHR_Rthm,'LineWidth',1.5,'Color',"#386cb0",'DisplayName','Valid Rhythm');
    hold on
    ylim([0 1.25]) % needed different scaling of y axis due to 0s and 1s from invalids
    xlim([CTDmin CTDmax]) % empty gap from 200 - 300ms doesn't look good w/ automatic limits; needs to be hardcoded since respective var's are inside loop
    yticks(0:0.1:1)
    xlabel('Cue-Target Delay (ms)')
    ylabel('Hit Rate')
    title(dateStr,' Performance as function of cue-target delay')
    legend(line,'Location','northwest','FontSize',8)
    legend('boxoff')
    txt         = string('number of trials: '+string(NumTrls));
    text(1200,1.10,txt,"HorizontalAlignment","right",FontWeight='bold')
    hold off;
    %TILE 2
    nexttile
    
    line(1)                 = plot(Bins,InvDO_HR_RthmDtnd,'LineWidth',1.1,'Color',"#fdc086",'DisplayName','Invalid DO Rhythm');
    hold on;
    line(2)                 = plot(Bins,InvSO_HR_RthmDtnd,'LineWidth',1.1,'Color',"#7fc97f",'DisplayName','Invalid SO Rhythm');
    hold on;
    line(3)                 = plot(Bins,ValHR_RthmDtnd,'LineWidth',1.5,'Color',"#386cb0",'DisplayName','Valid Rhythm');
    ylim([-1 1])
    xlim([CTDmin CTDmax]) % empty gap from 200 - 300ms doesn't look good w/ automatic limits; needs to be hardcoded since respective var's are inside loop
    xlabel('Cue-Target Delay (ms)')
    ylabel('Percentage above/below mean H.R.')
    title(string(dateStr),'Detrended performance as function of cue-target delay')
    legend(line,'Location','northwest','FontSize',8)
    legend('boxoff')
    text(1200,0.8,txt,"HorizontalAlignment",'right',FontWeight='bold')
    %tile3
    nexttile
    Fs          = 1000/slide;%fft function is based off milliseconds, highest SAMPLING FREQUENCY we can get though is 100 due to the binning of the slide (10 ms)
    T           = 1/Fs; %SAMPLING PERIOD
    L           = length(CtrTmeWindow); 
    timeVec     = (0:L-1)*T; % TIME VECTOR
    RthmTrnsfrm = fft(ValHR_RthmDtnd,L);
    P2          = abs(RthmTrnsfrm/L);
    P1          = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f           = Fs*(0:(L/2))/L;
    plot(f,P1,'Color',"#386cb0")
    title("Amplitude spectrum of HitRate(t)")
    xlabel("f (Hz)")
    ylabel("Abs. Proportion of freq.")
    title(string(dateStr)," Amplitude spectrum of HitRate(t)")    
    
    fileName    = string(dailyPlots+"\"+dateStr+"RhythmBehavior.png");
    fig         = gcf;
    saveas(fig,fileName);
    close(fig);
end
%8/18/2023 - ANYTHING BELOW THIS LINE IS GOING TO BROUGHT UP OR REPLACED AS I COMPLETE THE ACTUAL FUNCTION (ANYTHING ABOVE THIS LINE)
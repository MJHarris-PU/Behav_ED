function [cumulative] = cumulative(Data)
%  CUMULATIVE FUNCTION ANALYZING BEHAVIORAL RHYTHMS OF PERFORMANCE - CUMULATIVE
% SOME GLOBAL VARIABLES
minRT           = 150;
maxRT           = 800;
%PABLO FILE LOC'S
%cumulativePlots = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\CUMULATIVE';
%cumulativePlots = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\CUMULATIVE';
%FREDDIE FILE LOC'S
cumulativePlots = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\FreddieCUMULATIVE';
%cumulativePlots = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS\FIGURES\FreddieCUMULATIVE';
%CALCULATE THE PERFORMANCE FOR A GIVEN SESSION
sessions        = unique(Data.date);
CTDmin          = 500; % CHANGES DEPENDING ON CONDITION
CTDmax          = 1400;
CTDstep         = 10;
Binwidth        = 50;
Bins            = ((CTDmin+(CTDmin+Binwidth))/2):CTDstep:((CTDmax+(CTDmax-Binwidth))/2);
NumBins         = length(Bins);
ValHRarray      = zeros(height(sessions),NumBins);
InvSOHRarray    = zeros(height(sessions),NumBins);
InvDOHRarray    = zeros(height(sessions),NumBins);
ValHRarrayTnd   = zeros(height(sessions),NumBins);
InvSOHRarrayTnd = zeros(height(sessions),NumBins);
InvDOHRarrayTnd = zeros(height(sessions),NumBins);
%zeros(1,length(dummyEndTme:dummySlide:dummyMaxTme));
NumTrls     = 0;
StrtDate    = datetime(sessions(1));
EndDate     = datetime(sessions(end));
StrtDateStr = string(StrtDate,'MM-dd-yyyy');
EndDateStr  = string(EndDate,'MM-dd-yyyy');
for i = 1:length(sessions)
    DailyTbl                = Data(Data.date == sessions(i),:);                 % TABLE TO CAPTURE THE CURRENT ITERATION (SESSION) 
    EyeBrkIdx               = ismember(DailyTbl.break_eye,'true');              % REMOVING EYE BREAKS FROM THE START
    PreCueLevBrkIdx         = ismember(DailyTbl.break_eye,'false') & ismember(DailyTbl.break_lever,'true') & (DailyTbl.t_cue_on == 0);  % LEVER BREAK TRIALS BEFORE A CUE ONSET; INDEPENDENT OF EyeBrkIdx 
    RmvIdx                  = logical(EyeBrkIdx + PreCueLevBrkIdx);
    UseableTbl              = DailyTbl;                                         % SEPARATE TABLE INSTEAD OF MODIFYING DailyTbl. I WANT TO KEEP DAILY TABLE INTACT, SO THAT I CAN TRACK THOSE REMOVE TRIALS LATER ON IF NEEDED
    UseableTbl(RmvIdx,:)    = [];
    NumTrls                 = NumTrls + height(UseableTbl);
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
    %ValHR_Rthm(isnan(ValHR_Rthm))       = 0;
    %InvSO_HR_Rthm(isnan(InvSO_HR_Rthm)) = 0;
    %InvDO_HR_Rthm(isnan(InvDO_HR_Rthm)) = 0;
    ValHR_RthmDtnd          = detrend(ValHR_Rthm,"omitnan");
    InvSO_HR_RthmDtnd       = detrend(InvSO_HR_Rthm,"omitnan");
    InvDO_HR_RthmDtnd       = detrend(InvDO_HR_Rthm,"omitnan");
    
    ValHRarrayTnd(i,:)      = ValHR_Rthm;           % ARRAY OF VALID HIT RATES FOR EACH SESSION PRE-DETREND
    InvSOHRarrayTnd(i,:)    = InvSO_HR_Rthm;        % ARRAY OF INV S.O. HR PER SESSION PRE-DETREND
    InvDOHRarrayTnd(i,:)    = InvDO_HR_Rthm;        % ARRAY OF INV D.O. HR PER SESSION PRE-DETREND
    
    ValHRarray(i,:)         = ValHR_RthmDtnd;       % ARRAY OF DETRENDED VALID HIT RATES FOR EACH SESSION (EACH ROW)
    InvSOHRarray(i,:)       = InvSO_HR_RthmDtnd;    % ARRAY OF DETRENDED INVALID SO HIT RATES FOR EACH SESSION (EACH ROW) 
    InvDOHRarray(i,:)       = InvDO_HR_RthmDtnd;    % ARRAY OF DETRENDED INVALID DO HIT RATES FOR EACH SESSION (EACH ROW)
end
% MEAN & SEM VALUES FOR EACH CONDITION OF H.R. PRE-DETREND
MnValHRtrnd         = mean(ValHRarrayTnd,"omitnan");
MnInvSOHRtrnd       = mean(InvSOHRarrayTnd,"omitnan");
MnInvDOHRtrnd       = mean(InvDOHRarrayTnd,"omitnan");
MnValHRtrndSEM      = std(ValHRarrayTnd,"omitnan")/sqrt(height(sessions));
MnInvSOHRtrndSEM    = std(InvSOHRarrayTnd,"omitnan")/sqrt(height(sessions));
MnInvDOHRtrndSEM    = std(InvDOHRarrayTnd,"omitnan")/sqrt(height(sessions));
% MEAN & SEM VALUES FOR EACH CONDITION OF HIT RATE DETRENDED
MnValHRarray        = mean(ValHRarray,"omitnan");
MnInvSOHRarray      = mean(InvSOHRarray,"omitnan");
MnInvDOHRarray      = mean(InvDOHRarray,"omitnan");
ValHRarraySEM       = std(ValHRarray,"omitnan")/sqrt(height(sessions));
InvSOHRarraySEM     = std(InvSOHRarray,"omitnan")/sqrt(height(sessions));
InvDOHRarraySEM     = std(InvDOHRarray,"omitnan")/sqrt(height(sessions));
% STRING INDICATING PERIOD OF VISUALIZATION
datestr     = strcat(StrtDateStr," thru ",EndDateStr);
    
%SUBPLOTS TILED LAYOUT 
%REPEATED VARIABLES
figure('Position', [1,1,2100,800]);
%t               = tiledlayout(1,3);
tiledlayout(1,3);
x2              = [Bins,fliplr(Bins)]; % array that holds x positions. Is used repeatedly below
ValColor        = [195 210 231]/256;
InvDOcolor      = [253 217 182]/256;
InvSOcolor      = [203 233 203]/256;
txt             = string('Total trials: '+string(NumTrls));
%TILE 1 - TRENDED GRAPH 
nexttile
inBetween       = [(MnInvDOHRtrnd+MnInvDOHRtrndSEM),fliplr(MnInvDOHRtrnd-MnInvDOHRtrndSEM)];
fill1           = fill(x2, inBetween, InvDOcolor, 'LineStyle','none');
fill1.FaceAlpha = 0.3;
hold on
line(1)         = plot(Bins,MnInvDOHRtrnd,'LineWidth',1.3,'Color','#fdc086','DisplayName','Invalid DO Rhythm');
hold on
inBetween       = [(MnInvSOHRtrnd+MnInvSOHRtrndSEM),fliplr(MnInvSOHRtrnd-MnInvSOHRtrndSEM)];
fill2           = fill(x2, inBetween, InvSOcolor, 'LineStyle','none');
fill2.FaceAlpha = 0.3;
hold on 
line(2)         = plot(Bins,MnInvSOHRtrnd,'LineWidth',1.3,'Color','#7fc97f','DisplayName','Invalid SO Rhythm');
hold on
inBetween       = [(MnValHRtrnd+MnValHRtrndSEM),fliplr(MnValHRtrnd-MnValHRtrndSEM)];
fill3           = fill(x2, inBetween, ValColor, 'LineStyle','none');
fill3.FaceAlpha = 0.7;
hold on
line(3)         = plot(Bins,MnValHRtrnd,'LineWidth',1.8,'Color','#386cb0','DisplayName','Valid Rhythm');
ylim([0 1.25])
xlim([CTDmin CTDmax])
yticks(0:0.1:1)
xlabel('Cue-Target Delay (ms)')
ylabel('Hit Rate')
title('Cumulative performance as function of CTD',datestr)
legend(line,'Location','northwest','FontSize',8)
legend('boxoff')
text(1200,1.1,txt,"HorizontalAlignment","right",FontWeight='bold')
hold off
%TILE 2 - DETRENDED GRAPH
nexttile
inBetween       = [(MnInvDOHRarray+InvDOHRarraySEM),fliplr(MnInvDOHRarray-InvDOHRarraySEM)];
fill1 = fill(x2, inBetween, InvDOcolor, 'LineStyle','none');
fill1.FaceAlpha = 0.3;
hold on;
line(1)         = plot(Bins, MnInvDOHRarray,'LineWidth',1.3,'Color','#fdc086','DisplayName','Invalid DO Rhythm');
hold on;
inBetween       = [(MnInvSOHRarray+InvSOHRarraySEM),fliplr(MnInvSOHRarray-InvSOHRarraySEM)];
fill2 = fill(x2, inBetween, InvSOcolor, 'LineStyle','none');
fill2.FaceAlpha = 0.3;
hold on;
line(2)         = plot(Bins, MnInvSOHRarray,'LineWidth',1.3,'Color','#7fc97f','DisplayName','Invalid SO Rhythm');
hold on;        
inBetween       = [(MnValHRarray+ValHRarraySEM),fliplr(MnValHRarray-ValHRarraySEM)];
fill3 = fill(x2, inBetween, ValColor, 'LineStyle','none');
fill3.FaceAlpha = 0.7;
hold on;
line(3)         = plot(Bins, MnValHRarray, 'LineWidth',1.8,'Color',"#386cb0",'DisplayName','Valid Rhythm');
ylim([-0.2 0.2])
xlim([CTDmin CTDmax])
yticks(-0.2:0.05:0.2)
xlabel('Cue-target delay (ms)')
ylabel('Hit Rate % above/below mean')
title('Cumulative detrended performance as function of CTD',datestr)
legend(line,'Location','northwest','FontSize',8)
legend('boxoff')
text(1200,0.15,txt,"HorizontalAlignment","right",FontWeight='bold')
hold off
% TILE 3 - FREQUENCY SPECTRUM DISTRIBUTION
nexttile
Fs              = 1000/CTDstep;
%T               = 1/Fs;
L               = length(Bins);
%t               = (0:L-1)*T;
RthmTrnsfrm     = fft(MnValHRarray,L);
%HannL           = hann(L,"periodic");       %COMMENT OUT IF NOT WORKING AFTER A WHILE
%RthmTrnsfrm     = fft(hann(L,"periodic").*MnValHRarray,L);  %COMMENT OUT IF NOT WORKING AFTER A WHILE 
% ^^^ NOT GETTING THE HANNING WINDOW YET ^^^
P2              = abs(RthmTrnsfrm/L);
P1              = P2(1:L/2+1);
P1(2:end-1)     = 2*P1(2:end-1);     
f               = Fs*(0:(L/2))/L;
plot(f,P1,'Color',"#386cb0")
title("Amplitude spectrum of HitRate(t)")
xlabel("f (Hz)")
ylabel("Abs. Proportion of freq.")
title("Amplitude spectrum of HitRate(t) ",datestr)
% SAVE
fileName    = string(cumulativePlots+"\"+"CumulativeRhythm.png");
cumulative  = gcf;
saveas(cumulative,fileName);
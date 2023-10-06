function [PerformanceMetrics] = PerfMetric(Data)
%  PERFMETRIC FUNCTION THAT CALCULATES METRICS ON ANY SUBTABLE OF 'DATA'
PerformanceMetrics = table();
minRT = 150;
maxRT = 800;
% DAILY OVERALL PERFORMANCE
sessions = unique(Data.date);
    for i = 1:length(sessions)
        PerformanceMetrics.Date(i) = sessions(i);
        
        DailyTbl         = Data(Data.date == sessions(i),:); %Table for particular day in 'sessions'
        
        %SECTION ONE - REMOVE TRIALS NOT ANALYZED
        % Index for eye breaks and instances where cue never comes up
        % before lever breaks. These should not be counted towards
        % trial counts...
        EyeBrkIndx          = ismember(DailyTbl.break_eye, 'true'); %Eye breaks should be removed
        PreCueLevBrksIndx   = ismember(DailyTbl.break_eye,'false') & ismember(DailyTbl.break_lever,'true') & (DailyTbl.t_cue_on==0); % Break levers before cue onset. Needs to be exlcusive of EyeBrks so both logical arrays can be used together below.
        
        %RELATED COUNTS
        NumEyeBrk           = sum(EyeBrkIndx); %Number of eye breaks
        NumPreCueLevBrk     = sum(PreCueLevBrksIndx); %Number of lever break before cue onset
        NumTrls             = height(DailyTbl)-(sum(EyeBrkIndx+PreCueLevBrksIndx)); %Number of trials once removing above exceptions
        %BEGIN REMOVING THE UNUSEABLE TRIALS 
        RmvIdx              = logical(EyeBrkIndx + PreCueLevBrksIndx); % Makes index of trials to remove
        NumRmvTrls          = sum(RmvIdx); %Number of removable trials
        
        %NEW TABLE EXCLUSING THE REMOVEABLE TRIALS
        UseableTbl          = DailyTbl; %Copy table initially
        UseableTbl(RmvIdx,:)= []; %Change the secondary table to actually remove unuseable trials
        
        %SECTION TWO - MAKE VARIOUS INDICES THAT CALCULATE THE BASIC METRICS
        %(Overall Trials, Correct, Valid, Invalid, Catch, Early Break, Eye Break, Lever Break, etc.) 
        %Start w/ valids, invalids, catches, then base other measures off
        %aggregates of these
    
        ValIdx              = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 1:8; %Daily valid trial index
        ValIdx              = logical(sum(ValIdx,2)); %Transforms the index from 8 columns to one; should only have "1" for each row
        ValTrlTbl           = UseableTbl(ValIdx,:); %Table for valid trials
        ValCorr             = ((ValTrlTbl.t_contrast_change + minRT) <= ValTrlTbl.t_release) & (ValTrlTbl.t_release <= (ValTrlTbl.t_contrast_change + maxRT)); %Index for correct valid trials
        ValIncorr           = ValTrlTbl(ValCorr == 0, :); % Table showing the incorrect valid trials
        NumValCorr          = sum(ValCorr); %Number of valid correct trials
        ValCorrTbl          = ValTrlTbl(ValCorr,:); %Table for valid trial performance
        %Invalid and Catch variables receive similar descriptions as valids
        InvalSOIdx          = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 9;
        InvalSOTrlTbl       = UseableTbl(InvalSOIdx,:);
        InvalSOCorr           = ((InvalSOTrlTbl.t_contrast_change + minRT) <= InvalSOTrlTbl.t_release) & (InvalSOTrlTbl.t_release <= InvalSOTrlTbl.t_contrast_change + maxRT);
        InvalSOIncorr         = InvalSOTrlTbl(InvalSOCorr == 0,:);
        NumInvalSOCorr        = sum(InvalSOCorr);
        InvalSOCorrTbl        = InvalSOTrlTbl(InvalSOCorr,:);
        InvalDOIdx            = UseableTbl.if_catch_trial == 2 & UseableTbl.if_valid_trial == 10; 
        InvalDOTrlTbl         = UseableTbl(InvalDOIdx,:);
        InvalDOCorr           = ((InvalDOTrlTbl.t_contrast_change + minRT) <= InvalDOTrlTbl.t_release) & (InvalDOTrlTbl.t_release <= InvalDOTrlTbl.t_contrast_change + maxRT);
        InvalDOIncorr         = InvalDOTrlTbl(InvalDOCorr == 0,:);
        NumInvalDOCorr        = sum(InvalDOCorr);
        InvalDOCorrTbl        = InvalDOTrlTbl(InvalDOCorr,:);
        CatchIdx            = UseableTbl.if_catch_trial == 1;
        CatchTrlTbl         = UseableTbl(CatchIdx,:);
        CatchCorr           = ismember(CatchTrlTbl.break_lever,'false');
        CatchCorrTbl        = CatchTrlTbl(CatchCorr,:);
        CatchIncorr         = ismember(CatchTrlTbl.break_lever,'true');
        CatchIncorrTbl      = CatchTrlTbl(CatchIncorr,:);
        NumCatchCorr        = sum(CatchCorr);
        NumCatchIncorr      = sum(CatchIncorr);
        %ADDITIONAL METRICS FROM THE ABOVE TABLES/VARS TO USE
        NumCorr             = sum(ValCorr) + sum(InvalSOCorr) + sum(InvalDOCorr) + sum(CatchCorr); %Number of corrects
        EarlyLateResp       = [ValIncorr(ismember(ValIncorr.break_lever,'false'),:); InvalSOIncorr(ismember(InvalSOIncorr.break_lever,'false'),:);InvalDOIncorr(ismember(InvalDOIncorr.break_lever,'false'),:)]; %Trials that arent lever breaks but were incorrect (out of response windows)
        LevBrkTrls          = [ValIncorr(ismember(ValIncorr.break_lever,'true'),:); InvalSOIncorr(ismember(InvalSOIncorr.break_lever,'true'),:); InvalDOIncorr(ismember(InvalDOIncorr.break_lever,'true'),:) ; CatchTrlTbl(ismember(CatchTrlTbl.break_lever,'true'),:)]; %Trials that are part of the wrong category
        PropCorr            = NumCorr/NumTrls; %not counting eye breaks or breaks before cue appears toward correct proportion (for now. follow back up on this) 
        %OTHER ASPECTS OF PERFORMANCE (RT, etc...)
        ValTrlMnRT          = mean(ValCorrTbl.t_release - ValCorrTbl.t_contrast_change); % Mean RT for correct valid trials
        InvalSOTrlMnRT        = mean(InvalSOCorrTbl.t_release - InvalSOCorrTbl.t_contrast_change); % Mean RT for correct invalid SO trials
        InvalDOTrlMnRT        = mean(InvalDOCorrTbl.t_release - InvalDOCorrTbl.t_contrast_change); % Mean RT for correct invalid DO trials
        PropValCorr         = NumValCorr/height(ValTrlTbl); % Proportion of Valid Corrects
        PropInvalSOCorr       = NumInvalSOCorr/height(InvalSOTrlTbl);
        PropInvalDOCorr       = NumInvalDOCorr/height(InvalDOTrlTbl); % Proportion of correct invalids
        PropCatchCorr       = NumCatchCorr/(NumCatchCorr+NumCatchIncorr); % Proportion of correct catches
        % BEGIN FILLING IN PERFORMANCE METRICS TABLE             
        PerformanceMetrics.TotalTrialCount(i)       = NumTrls;
        PerformanceMetrics.Correct(i)               = sum(NumCorr);
        PerformanceMetrics.ValidTrials(i)           = height(ValTrlTbl);
        PerformanceMetrics.InvalidSOTrials(i)         = height(InvalSOTrlTbl);
        PerformanceMetrics.InvalidDOTrials(i)         = height(InvalDOTrlTbl);
        PerformanceMetrics.CatchTrials(i)           = height(CatchTrlTbl);
        PerformanceMetrics.EyeBreaks(i)             = sum(NumEyeBrk); 
        PerformanceMetrics.EarlyLateResponse(i)     = height(EarlyLateResp);
        PerformanceMetrics.LeverBreaks(i)           = height(LevBrkTrls);
        PerformanceMetrics.OverallPropCorrect(i)    = PropCorr;
        PerformanceMetrics.ValidPropCorrect(i)      = PropValCorr;
        PerformanceMetrics.InvalidSOPropCorrect(i)    = PropInvalSOCorr; 
        PerformanceMetrics.InvalidDOPropCorrect(i)    = PropInvalDOCorr;
        PerformanceMetrics.CatchPropCorrect(i)      = PropCatchCorr;
        PerformanceMetrics.ValidRT(i)               = ValTrlMnRT;
        PerformanceMetrics.InvalidSORT(i)             = InvalSOTrlMnRT;
        PerformanceMetrics.InvalidDORT(i)             = InvalDOTrlMnRT;
                
        % QUADRANT BASED METRICS - FOR LOOP THAT RUNS THROUGH THE DIFFERENT
        % QUADRANTS & GETS THE REGIONS-SPECIFIC PERFORMANCE VALUES
        Quads               = {[-1 1],[1 1],[-1 -1],[1 -1]}; % goes in this order: UL, UR, LL, LR; first value cue_xloc, second value cue_yloc; when calling to specific quad, use syntax 'quads{1}' & 'quads{1}(1)'
        QuadLabels          = ["UL","UR","LL","LR"]; % VERY IMPORTANT THAT THESE ALIGN WITH THE ORDER OF THE X & Y's IN Quads ...
       
        for j = 1:length(Quads)
            
            QuadValTbl          = ValTrlTbl(ValTrlTbl.cue_x_loc == Quads{j}(1) & ValTrlTbl.cue_y_loc == Quads{j}(2),:);
            QuadCorrValTbl      = ValCorrTbl(ValCorrTbl.cue_x_loc == Quads{j}(1) & ValCorrTbl.cue_y_loc == Quads{j}(2),:);
            QuadInvalSOTbl        = InvalSOTrlTbl(InvalSOTrlTbl.cue_x_loc == Quads{j}(1) & InvalSOTrlTbl.cue_y_loc == Quads{j}(2),:);
            QuadCorrInvalSOTbl    = InvalSOCorrTbl(InvalSOCorrTbl.cue_x_loc == Quads{j}(1) & InvalSOCorrTbl.cue_y_loc == Quads{j}(2),:);
            
            QuadInvalDOTbl        = InvalDOTrlTbl(InvalDOTrlTbl.cue_x_loc == Quads{j}(1) & InvalDOTrlTbl.cue_y_loc == Quads{j}(2),:);
            QuadCorrInvalDOTbl    = InvalDOCorrTbl(InvalDOCorrTbl.cue_x_loc == Quads{j}(1) & InvalDOCorrTbl.cue_y_loc == Quads{j}(2),:);
            QuadCatchTbl        = CatchTrlTbl(CatchTrlTbl.cue_x_loc == Quads{j}(1) & CatchTrlTbl.cue_y_loc == Quads{j}(2),:);    
            QuadCorrCatchTbl    = CatchCorrTbl(CatchCorrTbl.cue_x_loc == Quads{j}(1) & CatchCorrTbl.cue_y_loc == Quads{j}(2),:);
            Quad_ValCorr        = append(QuadLabels(j),"_ValidCorrect");
            Quad_InvalSOCorr      = append(QuadLabels(j),"_Invalid_SO_Correct");
            Quad_InvalDOCorr      = append(QuadLabels(j),"_Invalid_DO_Correct");
            Quad_CatchCorr      = append(QuadLabels(j),"_CatchCorrect");
            Quad_ValMnRT        = append(QuadLabels(j),"_ValidMeanRT");
            Quad_InvalSOMnRT      = append(QuadLabels(j),"_Invalid_SO_MeanRT");
            Quad_InvalDOMnRT      = append(QuadLabels(j),"_Invalid_DO_MeanRT");
            PerformanceMetrics.(Quad_ValCorr)(i)    = height(QuadCorrValTbl)/height(QuadValTbl);
            PerformanceMetrics.(Quad_InvalSOCorr)(i)  = height(QuadCorrInvalSOTbl)/height(QuadInvalSOTbl);
            PerformanceMetrics.(Quad_InvalDOCorr)(i)  = height(QuadCorrInvalDOTbl)/height(QuadInvalDOTbl);
            PerformanceMetrics.(Quad_CatchCorr)(i)  = height(QuadCorrCatchTbl)/height(QuadCatchTbl);
            PerformanceMetrics.(Quad_ValMnRT)(i)    = mean(QuadCorrValTbl.t_release - QuadCorrValTbl.t_contrast_change);
            PerformanceMetrics.(Quad_InvalSOMnRT)(i)  = mean(QuadCorrInvalSOTbl.t_release - QuadCorrInvalSOTbl.t_contrast_change);
            PerformanceMetrics.(Quad_InvalDOMnRT)(i)  = mean(QuadCorrInvalDOTbl.t_release - QuadCorrInvalDOTbl.t_contrast_change);
        end   
    end
% barY = [mean(PerformanceMetrics.ValidPropCorrect,"omitnan"),mean(PerformanceMetrics.InvalidSOPropCorrect,"omitnan"),mean(PerformanceMetrics.InvalidDOPropCorrect,"omitnan"),mean(PerformanceMetrics.CatchPropCorrect,"omitnan")];
% ValColor        = [195 210 231]/256;
% InvDOcolor      = [253 217 182]/256;
% InvSOcolor      = [203 233 203]/256;
% Catchcolor      = [255,255,165]/256;
% 
% PerfBar = bar(barY);
% title("Hit Rates Per Trial Condition")
% xticklabels({"Valid","Invalid S.O.","Invalid D.O.","Catch"});
% PerfBar.FaceColor = 'flat';
% PerfBar.CData(1,:) = ValColor;
% PerfBar.CData(2,:) = InvSOcolor;
% PerfBar.CData(3,:) = InvDOcolor;
% PerfBar.CData(4,:) = Catchcolor;
%
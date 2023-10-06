function [Object] = Object(Data) % I think the return value will be a structure, containing two tables: horizontal and vertical
%  OBJECT FUNCTION THE COMPARES SAME & DIFFERENT OBJECT PERFORMANCE
Object          = struct();
bar_ori         = [1,2]; % 1 == Vertical; 2 == Horizontal
minRT           = 150;
maxRT           = 800;
sessions = unique(Data.date);
    for i = 1:length(sessions)
        ObjTbl.Date(i)      = sessions(i);
        DailyTbl            = Data(Data.date == sessions(i),:); %Table for particular day in 'sessions'
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
% DETERMINED THAT VRT AND HRZ BAR SEPARATION ISN'T FULLY NECESSARY FOR OBJECT BASED, HOWEVER IT MAY COME IN USEFUL, JUST NEED TO SEPARATE ALONG EITHER 'if_valid_trial' == 9 OR 10
        %USING THE USEABLE TABLE, SEPARATE TRIALS INTO HORIZONTAL AND VERTICAL TABLES
        VrtIdx              = UseableTbl.bar_ori == 1;
        VrtTbl              = UseableTbl(VrtIdx,:);
        HrzIdx              = UseableTbl.bar_ori == 2;
        HrzTbl              = UseableTbl(HrzIdx,:);
 
        %ADDITIONAL SUBTABLES FOR EACH VRT & HRZ TO SEPARATE VALID, INVALID S.O., & INVALID D.O. CONDITIONS
        VrtValTblIdx        = VrtTbl.if_catch_trial == 2 & VrtTbl.if_valid_trial == 1:8;
        VrtValTbl           = VrtTbl(logical(sum(VrtValTblIdx,2)),:);
        VrtInvalTblSO       = VrtTbl((VrtTbl.if_catch_trial == 2 & VrtTbl.if_valid_trial == 9),:);
        VrtInvalTblDO       = VrtTbl((VrtTbl.if_catch_trial == 2 & VrtTbl.if_valid_trial == 10),:);
        HrzValTblIdx        = HrzTbl.if_catch_trial == 2 & HrzTbl.if_valid_trial == 1:8;
        HrzValTbl           = HrzTbl(logical(sum(HrzValTblIdx,2)),:);
        HrzInvalTblSO       = HrzTbl((HrzTbl.if_catch_trial == 2 & HrzTbl.if_valid_trial == 9),:);
        HrzInvalTblDO       = HrzTbl((HrzTbl.if_catch_trial == 2 & HrzTbl.if_valid_trial == 10),:);
% 8/10/23 POSTPONING CONTINUING OBJECT FOR NOW AND SWITCHING TO RHYTHM
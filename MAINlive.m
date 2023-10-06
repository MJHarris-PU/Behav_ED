%% 
% FIRST VERSION OF ANALYSIS CODE FOR MPH EXPERIMENT

% Setting up directories for data file import and useful libraries
clear;
clc;
% Home Desktop - ONLY ONE USED AT A TIME
%dataDir = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\DATA';
%homeDir = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS';


%Work Laptop -ONLY ONE USED AT A TIME
dataDir = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\DATA\';

% CURRENTLY ASSIGNED TO FREDDIES DATA FOLDER
%dataDir = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\DATAF\';
%dataDir = 'C:\Users\mjh\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\DATAF\';



homeDir = 'C:\Users\mjhar\OneDrive - Princeton University\PU_RESEARCH_SPECIALIST\!KASTNER_LAB\PHARMACOLOGICAL_MANIPULATION\ANALYSIS';
addpath(homeDir);
% .edf File Converter Location
%% 
% READ IN DATA FILES AND ARRANGE, CONCATENATE TO MAIN DATA TABLE, ORGANIZE AS 
% NECESSARY

cd(dataDir); % Enter Data Folder
DataFldr = dir("*.csv"); % Capture File Details
Files = {DataFldr.name}; % Cell Array of File Names
Data  = cell(1,numel(DataFldr)); % Pre-allocated space for table equal to size of data files

% FOR LOOP THROUGH FILES, LOADS CSV INTO TABLE, MAKE VARIABLE OF FIRST SIX INTS OF
% CURRENT FILE NAME, ADDS CSV TO TABLE + COLUMN WITH THE DATE
for i=1:length(Files)
    Data{i} = readtable(Files{i}); 
    tmpDate = string(Files{i}(1:6));
    Data{i}.date(:) = tmpDate;
end

Data = cat(1,Data{:}); %Concatenate cell array of matrices into large table
Data = [Data(:,end) Data(:,1:end-1)]; %reorder so date is first column
Data.date = datetime(Data.date,'InputFormat',"yyMMdd"); %change format of date column to datetime
%% 
% DEFINE SUBCATEGORIES (SUBTABLES) THAT ARE NECESSARY FOR LATER OBJ-BASED ANALYSIS 

% % Catch Trial Table
% CatchIdx = Data.if_catch_trial == 1; %indices that are catches
% Catch = Data(CatchIdx,:); %New table according to this index
% 
% % Invalid Table
% InvIdx = Data.if_catch_trial == 2 & Data.if_valid_trial == 9 | Data.if_catch_trial == 2 & Data.if_valid_trial == 10; %If either 9 or 10 and not a catch
% Invalid = Data(InvIdx,:);
% 
% % Valid Table 
% vIdx = Data.if_catch_trial==2 & Data.if_valid_trial == 1:8; % If not a catch trial and not invalid
% vIdx = logical(sum(vIdx,2)); % sum across columns to make into single column. transform from double type to logical (T/F)
% Valid = Data(vIdx,:); %New table according to the generated index values
% 
%%
%TESTING/TROUBLESHOOTING SECTION
minimumRT = 150;
maximumRT = 800;

EyeBreaksIndex = ismember(Data.break_eye, 'true'); %find indices for eye breaks
PreCueLeverBreaksIndex   = ismember(Data.break_lever,'true') & ismember(Data.break_eye,'false') & (Data.t_cue_on==0);%find indices where lever broke before cue

NumberEyeBreaks = sum(EyeBreaksIndex);
NumberPreCueBreaks = sum(PreCueLeverBreaksIndex);

NumberTrials = height(Data)-(sum(EyeBreaksIndex+PreCueLeverBreaksIndex)); %get the number of trials that don't include the two above exceptions
RemoveIndex              = logical(EyeBreaksIndex + PreCueLeverBreaksIndex); %turn those trials into an index that can be removed
RemoveableTrials = sum(RemoveIndex);

UseableTable = Data; %make a secondary table that you will remove the unuseable trials from
UseableTable(RemoveIndex,:) = []; %change the secondary table to actually remove them

ValidIndex = UseableTable.if_catch_trial==2 & UseableTable.if_valid_trial == 1:8;
ValidIndex        = logical(sum(ValidIndex,2));
ValidTrialTable = UseableTable(ValidIndex,:);
ValidCorrects = ((ValidTrialTable.t_contrast_change + minimumRT) <= ValidTrialTable.t_release) & (ValidTrialTable.t_release <= (ValidTrialTable.t_contrast_change + maximumRT));
ValidIncorrects = ValidTrialTable(ValidCorrects==0,:);
NumberValidCorrects = sum(ValidCorrects);
ValidCorrectTable = ValidTrialTable(ValidCorrects,:);


InvalidIndex = (UseableTable.if_catch_trial == 2 & UseableTable.if_valid_trial == 9) | (UseableTable.if_catch_trial == 2 & UseableTable.if_valid_trial == 10);
InvalidTrialTable = UseableTable(InvalidIndex,:);
InvalidCorrects = ((InvalidTrialTable.t_contrast_change + minimumRT) <= InvalidTrialTable.t_release) & (InvalidTrialTable.t_release <= (InvalidTrialTable.t_contrast_change + maximumRT));
InvalidIncorrects = InvalidTrialTable(InvalidCorrects==0,:);
NumberInvalidCorrects = sum(InvalidCorrects);
InvalidCorrectTable = InvalidTrialTable(InvalidCorrects,:);

CatchIndex = UseableTable.if_catch_trial == 1;
CatchTrialTable = UseableTable(CatchIndex,:);
CorrectCatches = ismember(CatchTrialTable.break_lever,'false');
CorrectCatchTable = CatchTrialTable(CorrectCatches,:);
IncorrectCatches = ismember(CatchTrialTable.break_lever,'true');
IncorrectCatchTable = CatchTrialTable(IncorrectCatches,:);
NumberCorrectCatches = sum(CorrectCatches);
NumberIncorrectCatches = sum(IncorrectCatches);

%Valid & Invalid trials that are incorrect
EarlyLateResponses = [ValidIncorrects(ismember(ValidIncorrects.break_lever,'false'),:); InvalidIncorrects(ismember(InvalidIncorrects.break_lever,'false'),:)];


% 8/1/23: USE INCORRECT INDICES TO LOCATE THE BREAK LEVER TRIALS AND SUM
% ACROSS THE DIFFERENT CONDITIONS TO EQUAL THE LEVER BREAK TRIALS COLUMN
% MAY WANT TO SPLIT EARLY & LATE BREAK TRIALS? (BOTH STILL COMPRISE ERRORS)
% OR YOU HAVE AN OUTSIDE OF WINDOW INCORRECT, OR A BREAK LEVER, AND BOTH
% SHOULD EQUAL THE SUM OF INCORRECTS....SHOULD THEN MAKE AN INCORRECT
% COLUMN TOO...



% MAKE CORRECT BELOW AFTER VALID, INVALID, CATCH SEGMENTATIONS...
%ACTUAL Trial Indices:
% CorrectTrials = ((UseableTable.t_contrast_change + minimumRT) <= UseableTable.t_release) & (UseableTable.t_release <= (UseableTable.t_contrast_change + maximumRT)) | (UseableTable.if_catch_trial==1 & ismember(UseableTable.break_lever,'false'));

% JUST THE SUM OF VALUES BETWEEN DIFFERENT CONDITIONS (SHOULD EQUAL THE
% SAME IF SUMMING CORRECT TRIALS):
NumberCorrectTrials = sum(ValidCorrects) + sum(InvalidCorrects) + sum(CorrectCatches);

% Use the incorrect trials and catches that were lever breaks in a new
% table to capture all "lever breaks". Can then use the length of this as
% the value to input for lever break trials
LeverBreakTrials = [ValidIncorrects(ismember(ValidIncorrects.break_lever,'true'),:); InvalidIncorrects(ismember(InvalidIncorrects.break_lever,'true'),:); CatchTrialTable(ismember(CatchTrialTable.break_lever,'true'),:)];
NumberLeverBreakTrials = height(LeverBreakTrials);


IncorrectTrials = [ValidIncorrects;InvalidIncorrects;IncorrectCatchTable];
NumberIncorrectTrials = height(IncorrectTrials);
%% 
% 8/5/23 - MORE DESCRIPTIVE PERF METRICS

%PROP CORRECT
PropValCorr             = NumberValidCorrects/height(ValidTrialTable);
PropInvalCorr           = NumberInvalidCorrects/height(InvalidTrialTable);
PropCatchCorr           = NumberCorrectCatches/(NumberCorrectCatches+NumberIncorrectCatches);

%OVERALL RT VAL/INVAL
ValidTrialMeanRT        = mean(ValidCorrectTable.t_release - ValidCorrectTable.t_contrast_change);
InvalidTrialMeanRT      = mean(InvalidCorrectTable.t_release - InvalidCorrectTable.t_contrast_change);

%QUADRANT-BASED CORRECTS & RTs (OVERALL, UL, UR, LL, LR)
quads = {[-1 1],[1 1],[-1 -1],[1 -1]}; % goes in this order: UL, UR, LL, LR; first value is cue_xloc, second value is cue_yloc; when calling to specific quad, use syntax 'quads{1}' & 'quads{1}(1)'

%%
%FOR LOOP THAT WILL POPULATE THE QUADRANT LEVEL METRICS


QuadLabels = ["UL","UR","LL","LR"];
testTable = table();

for j = 1:length(quads)
    
    QuadValTbl          = ValidTrialTable(ValidTrialTable.cue_x_loc == quads{j}(1) & ValidTrialTable.cue_y_loc == quads{j}(2),:);
    QuadCorrectValTbl   = ValidCorrectTable(ValidCorrectTable.cue_x_loc == quads{j}(1) & ValidCorrectTable.cue_y_loc == quads{j}(2),:);
    QuadInvalTbl        = InvalidTrialTable(InvalidTrialTable.cue_x_loc == quads{j}(1) & InvalidTrialTable.cue_y_loc == quads{j}(2),:);
    QuadCorrectInvalTbl = InvalidCorrectTable(InvalidCorrectTable.cue_x_loc == quads{j}(1) & InvalidCorrectTable.cue_y_loc == quads{j}(2),:);
    QuadCatchTbl        = CatchTrialTable(CatchTrialTable.cue_x_loc == quads{j}(1) & CatchTrialTable.cue_y_loc == quads{j}(2),:);    
    QuadCorrectCatchTbl = CorrectCatchTable(CorrectCatchTable.cue_x_loc == quads{j}(1) & CorrectCatchTable.cue_y_loc == quads{j}(2),:);

    
    Quad_ValCorr        = append(QuadLabels(j),"_ValidCorrect");
    Quad_InvalCorr      = append(QuadLabels(j),"_InvalidCorrect");
    Quad_CatchCorr      = append(QuadLabels(j),"_CatchCorrect");
    Quad_ValMnRT        = append(QuadLabels(j),"_ValidMeanRT");
    Quad_InvalMnRT      = append(QuadLabels(j),"_InvalidMeanRT");

    
    testTable.(Quad_ValCorr)   = height(QuadCorrectValTbl)/height(QuadValTbl);
    testTable.(Quad_InvalCorr) = height(QuadCorrectInvalTbl)/height(QuadInvalTbl);
    testTable.(Quad_CatchCorr) = height(QuadCorrectCatchTbl)/height(QuadCatchTbl);
    testTable.(Quad_ValMnRT)   = mean(QuadCorrectValTbl.t_release - QuadCorrectValTbl.t_contrast_change);
    testTable.(Quad_InvalMnRT) = mean(QuadCorrectInvalTbl.t_release - QuadCorrectInvalTbl.t_contrast_change);
end
   
%% 
% SECTION FOR HORIZONTAL VS. VERTICAL BAR ORIENTATIONS. WILL BE USED IN OBJECT-BASED 
% ANALYSIS

VerticalIndex = UseableTable.bar_ori == 1;
VerticalTable = UseableTable(VerticalIndex,:);

VerticalValidTableIndex = VerticalTable.if_catch_trial == 2 & VerticalTable.if_valid_trial == 1:8;
VerticalValidTable      = VerticalTable(logical(sum(VerticalValidTableIndex,2)),:);
VerticalInvalidTableSO  = VerticalTable((VerticalTable.if_catch_trial == 2 & VerticalTable.if_valid_trial == 9),:);
VerticalInvalidTableDO  = VerticalTable((VerticalTable.if_catch_trial == 2 & VerticalTable.if_valid_trial == 10),:);

%% 
% SECTION FOR TESTING OUT RHYTHMIC RELATED MEASURES TO BE PLACED INTO RHYTHM 
% FUNCTION 

minimumRT = 150;
maximumRT = 800; 
EyeBreaksIndex              = ismember(Data.break_eye, 'true'); %find indices for eye breaks
PreCueLeverBreaksIndex      = ismember(Data.break_lever,'true') & ismember(Data.break_eye,'false') & (Data.t_cue_on==0);%find indices where lever broke before cue
RemoveIndex                 = logical(EyeBreaksIndex + PreCueLeverBreaksIndex);
UseableTable                = Data; %make a secondary table that you will remove the unuseable trials from
UseableTable(RemoveIndex,:) = []; %change the secondary table to actually remove them

% TRIALS THAT WERE BREAKLEVER == TRUE WILL NEED TO BE DROPPED SINCE THERE
% IS NO TIMING AVAILABLE TO CALCULATE THE CUE-TARGET INTERVAL




%VALID TABLE
ValidIndex                  = UseableTable.if_catch_trial == 2 & UseableTable.if_valid_trial == 1:8;
ValidIndex                  = logical(sum(ValidIndex,2));
ValidTable                  = UseableTable(ValidIndex,:);
ValidCorrects               = ((ValidTable.t_contrast_change + minimumRT) <= ValidTable.t_release) & (ValidTable.t_release <= (ValidTable.t_contrast_change + maximumRT));
ValidCorrectTable           = ValidTable(ValidCorrects,:);
ValidIncorrectTable         = ValidTable(ValidCorrects==0,:);
ValidIncorrectTable         = ValidIncorrectTable(ismember(ValidIncorrectTable.break_lever,'false'),:); % REMOVES LEVER BREAKS THAT CAN'T BE CALCULATED

%INVALID S.O. TABLE
InvalidSOIndex              = UseableTable.if_catch_trial == 2 & UseableTable.if_valid_trial == 9;
InvalidTableSO              = UseableTable(InvalidSOIndex,:);
InvalidSO_Corrects          = ((InvalidTableSO.t_contrast_change + minimumRT) <= InvalidTableSO.t_release) & (InvalidTableSO.t_release <= (InvalidTableSO.t_contrast_change + maximumRT));
InvalidSO_CorrectTable      = InvalidTableSO(InvalidSO_Corrects,:);
InvalidSO_IncorrectTable    = InvalidTableSO(InvalidSO_Corrects==0,:);
InvalidSO_IncorrectTable    = InvalidSO_IncorrectTable(ismember(InvalidSO_IncorrectTable.break_lever,'false'),:); % REMOVES LEVER BREAKS THAT CAN'T BE CALCULATED


%INVALID D.O. TABLE
InvalidDOIndex              = UseableTable.if_catch_trial == 2 & UseableTable.if_valid_trial == 10;
InvalidTableDO              = UseableTable(InvalidDOIndex,:);
InvalidDO_Corrects          = ((InvalidTableDO.t_contrast_change + minimumRT) <= InvalidTableDO.t_release) & (InvalidTableDO.t_release <= (InvalidTableDO.t_contrast_change + maximumRT));
InvalidDO_CorrectTable      = InvalidTableDO(InvalidDO_Corrects,:);
InvalidDO_IncorrectTable    = InvalidTableDO(InvalidDO_Corrects==0,:);
InvalidDO_IncorrectTable    = InvalidDO_IncorrectTable(ismember(InvalidDO_IncorrectTable.break_lever,'false'),:); % REMOVES LEVER BREAKS THAT CAN'T BE CALCULATED

%PICK UP HERE ^^^ NEED TO CHECK THESE. STILL NEED TO FILTER FOR CORRECT TRIALS (INCORRECTS WILL HAVE INACCURATE
%TIMES THAT CAN'T BE USED FOR THE TIME BINS..
Vcue_targ                  = ValidCorrectTable.t_contrast_change - ValidCorrectTable.t_cue_on;
VIcue_targ                 = ValidIncorrectTable.t_contrast_change - ValidIncorrectTable.t_cue_on;
ISOcue_targ                = InvalidSO_CorrectTable.t_contrast_change - InvalidSO_CorrectTable.t_cue_on;
ISOIcue_targ               = InvalidSO_IncorrectTable.t_contrast_change - InvalidSO_IncorrectTable.t_cue_on;
IDOcue_targ                = InvalidDO_CorrectTable.t_contrast_change - InvalidDO_CorrectTable.t_cue_on;
IDOIcue_targ               = InvalidDO_IncorrectTable.t_contrast_change - InvalidDO_IncorrectTable.t_cue_on;

%CAN ESTABLISH AN INITIAL RANGE, AND THEN FILTER TRIALS BASED ON WHAT
%CONTAINS A CUE-TARGET INTERVAL THAT IS WITHIN THE RANGE
startTime = 350; % CHANGE TO 350 ONCE CHECKED
endTime = 400; %CHANGE TO 400 ONCE CHECKED
slide   = 10;
maxTime = 1400;
%dur = (startTime:endTime);
%AS YOU GO THROUGH A FOR LOOP, FOR EACH FOLLOWING ITERATION, dur = dur + 10
%THIS BUMP THE RANGE UP, CHANGING THE TRIALS, AND RECALCULATING THE MEAN
%H.R.

% Perform for loop that creates indices based upon matching to the cue
% target intervals. Check the 
ValidHitRateRhythm              = [];
InvalSOHitRateRhythm            = [];
InvalDOHitRateRhythm            = [];
CenterTimeWindow                = [];
while endTime <= maxTime
    ValCorrects                 = sum(Vcue_targ >= startTime & Vcue_targ <= endTime);
    ValIncorrects               = sum(VIcue_targ >= startTime & VIcue_targ <= endTime);
    ValHitRate                  = ValCorrects/(ValCorrects+ValIncorrects);
    ValidHitRateRhythm(end+1)   = ValHitRate;
    
    InvalCorrectsSO             = sum(ISOcue_targ >= startTime & ISOcue_targ <= endTime);
    InvalIncorrectsSO           = sum(ISOIcue_targ >= startTime & ISOIcue_targ <= endTime);
    InvalSOHR                   = InvalCorrectsSO/(InvalCorrectsSO+InvalIncorrectsSO)
    InvalSOHitRateRhythm(end+1) = InvalSOHR;

    InvalCorrectsDO             = sum(IDOcue_targ >= startTime & IDOcue_targ <= endTime);
    InvalIncorrectsDO           = sum(IDOIcue_targ >= startTime & IDOIcue_targ <= endTime);
    InvalDOHR                   = InvalCorrectsDO/(InvalCorrectsDO+InvalIncorrectsDO);
    InvalDOHitRateRhythm(end+1) = InvalDOHR;
    
    CenterTimeWindow(end+1)     = (startTime+endTime)/2;
    startTime                   = startTime + slide;
    endTime                     = endTime+slide;
end

ValidHitRateRhythm(isnan(ValidHitRateRhythm))       = 0;
InvalSOHitRateRhythm(isnan(InvalSOHitRateRhythm))   = 0;
InvalDOHitRateRhythm(isnan(InvalDOHitRateRhythm))   = 0;

ValidHitRateRhythmDetrend       = detrend(ValidHitRateRhythm);
InvalidSOHitRateRhythmDetrend   = detrend(InvalSOHitRateRhythm);  
InvalidDOHitRateRhythmDetrend   = detrend(InvalDOHitRateRhythm);

%%
% PLOT FOR THE ABOVE
tmpDateTime = datetime(tmpDate,'InputFormat','yyMMdd');


line(1) = plot(CenterTimeWindow,ValidHitRateRhythm,'LineWidth',1.5,"DisplayName",'Valid Rhythm')
hold on
line(2) = plot(CenterTimeWindow,InvalDOHitRateRhythm,"DisplayName",'Invalid DO Rhythm')
hold on
line(3) = plot(CenterTimeWindow,InvalSOHitRateRhythm,"DisplayName",'Invalid SO Rhythm')
xlim([350 1400])
xlabel('Cue-Target Delay (ms)')
ylabel('Percentage above/below mean') % or just ('Hit Rate') if not detrended
title(string(tmpDateTime),'Performance as function of cue-target delay')
legend(line,'Location','southwest','FontSize',6)
legend('boxoff')
hold off
%Vcue_targ >= startTime & Vcue_targ <= endTime
%NOT DETRENDED YET, ROBER MENTIONED MAY HAVE BEEN AN OFF DAY WHICH IS WHY THE EXPECTED
%RAMPING CURVE IN PERFORMANCE IS NOT THERE...NEED TO LOOK FURTHER INTO THIS

%SIMPLE WAY TO PLOT THIS IS:
%'plot(CenterTimeWindow,ValidHitRateRhythmDetrend)'
%% 
% Fourier-related

Fs      = 100;      % sampling frequency (10 millisecond slide)
T       = 1/Fs;     % Sampling Period
%L       = length(CenterTimeWindow);  % Essentially "eighty-six one-hundreths of a second... I THINK THIS IS CORRECT...going with this since it's what allows me to run the calcs, signal length (end - start) = 860 ms divided by ten to go from 1 ms sample to 10 ms sample...
L       = length(CenterTimeWindow)*5
t       = (0:L-1)*T;%Time vector



%RhythmTransform     = fft(ValidHitRateRhythmDetrend);
RhythmTransform     = fft(ValidHitRateRhythmDetrend,L) %<< Smooths but changes frequency...I don't think correct
%RhythmTransform2     = [zeros(1,length(CenterTimeWindow)/2),RhythmTransform,zeros(1,length(CenterTimeWindow)/2)]; % also doesn't seem correct... 

P2 = abs(RhythmTransform/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1);%WILL NEED TO SEE IF THIS IS CORRECT. IT SORT OF LOOKS CORRECT, BUT...CAN'T SAY FOR SURE YET, GOING TO DO SOME MORE READING, ETC.
title("Amplitude spectrum of HitRate(t)")
xlabel("f (Hz)")
ylabel("Abs. Proportion of freq.")


%% 
% TESTING DIFFERENT COMBINATIONS OF ANY N# OF DAYS
% 
% 


%% 
% 7/27/23 PROBLEMS

currentTable = PerfMetric(Data)

CorrectIndex = ismember(Data.break_eye, 'false') & ((Data.t_contrast_change + 150) <= Data.t_release) & (Data.t_release <= (Data.t_contrast_change + 800));
CorrectTable = Data(CorrectIndex,:);

ValidIndex = Data.if_catch_trial==2 & Data.if_valid_trial == 1:8;
ValidIndex = logical(sum(ValidIndex,2));
ValidTable = Data(ValidIndex,:);

InvalidIndex = Data.if_catch_trial == 2 & Data.if_valid_trial == 9 | Data.if_catch_trial == 2 & Data.if_valid_trial == 10;
InvalidTable = Data(InvalidIndex,:);

CatchIndex = Data.if_catch_trial == 1;
CatchTable = Data(CatchIndex,:);

IncorrectIndex = find(CorrectIndex==0);
IncorrectTrials = Data(IncorrectIndex,:);

LeverBreakIndex = ismember(Data.break_lever,'true') & (Data.t_cue_on > 0) & (Data.t_release ~= Data.t_lever_break); %There were five additional trials here before adding Data.t_cue_on > 0 ... these would not represent a late release if the cue never even appeared.
LeverBreakTrials = Data(LeverBreakIndex,:); % count of these types was 76

LeverBreakCatchIdx = ismember(Data.break_lever,'true') & (Data.if_catch_trial == 1) & (Data.t_cue_on > 0) & (Data.t_release ~= Data.t_lever_break);
LeverBreakCatchTrials = Data(LeverBreakCatchIdx,:);

LeverBreakRegIdx = ismember(Data.break_lever,'true') & (Data.if_catch_trial == 2) & (Data.t_cue_on > 0) & (Data.t_release ~= Data.t_lever_break);
LeverBreakRegTrials = Data(LeverBreakRegIdx,:);

LateReleaseIndex = ismember(Data.break_eye,'false')&(Data.if_catch_trial==2)&(Data.t_release > (Data.t_contrast_change+800)); %<< include this towards performance calc.
LateReleaseTable = Data(LateReleaseIndex,:);

%DON'T COUNT TRIALS IN THE BELOW SECTION

EyeBreakIndex = ismember(Data.break_eye,'true');
EyeBreakTrials = Data(EyeBreakIndex,:);

PreCueLevBreakIndex = ismember(Data.break_lever,'true')&ismember(Data.break_eye,'false')& (Data.t_cue_on==0);
PreCueLevBreaks = Data(PreCueLevBreakIndex,:);


%EverythingElseIndex = 
%EverythingElseTable = 

%% 
% 


%% 
% STUFF CURRENTLY LISTED IN THIS SECTION MAY NOT BE RELEVANT WHEN USING FUNCTIONS

% % Quadrants << may not really be necessary when using function.
% quads = {[-1 1],[1 1],[-1 -1],[1 -1]}; %cell array that contains vectors of the four quadrants of E.D.
% [UL,UR,LL,LR] = quads{:}; % define each of the four quadrant variables
% % example line if I wanted to look at all valid cued trials that are upper
% % left quadrant:
% % Valid(Valid.cue_x_loc == UL(1) & Valid.cue_y_loc == UL(2),:) << **could
% % also work for tar_x_loc & tar_y_loc too...
%% 
% CAN LOOK AT ANALYSES HERE...EVENTUAL GOAL WILL BE TO MAKE SOME OF THESE IN 
% LINE SECTION THEIR OWN FUNCTIONS...

% % Proportion of correct trials (make applicable to whole data table or
% % subsets...)
% % ***Useful hint: for logicals expressed as 'true' 'false' in the break_eye
% % or break_lever columns, need to use ismember() function, instead of
% % logical operator...
% %minRT = 150;
% %rtLimit = 800; 
% 
% % STARTING W/ CORRECT
% %corrIdx = strcmp(Data.break_eye, Data.break_lever) % where both columns are false << using strcmp here wouldn't work if break_eye & break_lever could both potentially be true...FYI
% %corrIdx = strcmp(Data.break_eye, Data.break_lever) & Data.t_contrast_change + minRT < Data.t_release < Data.t_target_off + rtLimit; 
% % ^^^ CAN'T USE strcmp() since there are some break_eye and break_lever
% % columns that might both be 'true'
% %corrIdx = ismember(Valid.break_eye,'false') & ismember(Valid.break_lever, 'false') & Valid.t_contrast_change + minRT < Valid.t_release < Valid.t_target_off + rtLimit;
%% 
% AGGREGATED PERFORMANCE METRICS (ACROSS DAYS)
% 
% *** Currently something is quite messed up with how this calculates....need 
% to go back to Kacie script...

% minRT = 150; % minimum time following target onset to consider as a valid t_release time
% rtLimit = 800; % maximum time following end of target display to consider a valid release time
% 
% % VALID
% % Optional counts that affect proportion correct and should be removed
% % from denominator.
% EyeBrkIdx = ismember(Valid.break_eye,'true');
% EyeBrksum = nnz(EyeBrkIdx);
% 
% LtBrkIdx = ismember(Valid.break_eye,'false') & (Valid.t_release > (Valid.t_contrast_change + rtLimit));
% LtBrkSum = nnz(LtBrkIdx);
% 
% EarlyBrkIdx = ismember(Valid.break_eye,'false') & (Valid.t_release >= (Valid.t_contrast_change)) & (Valid.t_release <= (Valid.t_contrast_change + minRT));
% EarlyBrkSum = nnz(EarlyBrkIdx);
% 
% PreTargLevBrk = ismember(Valid.break_eye,'false') & (Valid.t_release >= (Valid.t_cue_on)) & (Valid.t_release < Valid.t_contrast_change);
% PreTargLevBrkSum = nnz(PreTargLevBrk);
% 
% PreCueLevBrk = ismember(Valid.break_eye, 'false') & (Valid.t_release < Valid.t_cue_on);
% PreCueLevBrkSum = nnz(PreCueLevBrk);
% 
% %6/8/23 PICK UP HERE ^^^ 
% 
% % Aggregate performance across all quadrants
% %corrIdx = ismember(Valid.break_eye,'false') & ismember(Valid.break_lever, 'false') & ((Valid.t_contrast_change + minRT) < Valid.t_release) & (Valid.t_release < (Valid.t_contrast_change + rtLimit));
% %corrIdx = ismember(Valid.break_eye,'false') & ((Valid.t_contrast_change + minRT) < Valid.t_release) & (Valid.t_release < (Valid.t_contrast_change + rtLimit)); %attempting w/o lever_break column...
% %corrIdx = ismember(Valid.break_eye,'false') & ((Valid.t_contrast_change + minRT) < Valid.t_lever_break) & (Valid.t_lever_break < (Valid.t_contrast_change + rtLimit)); %attempting w/o lever_break column...
% corrIdx = ismember(Valid.break_eye,'false') & ((Valid.t_contrast_change + minRT) < Valid.t_release) & (Valid.t_release < (Valid.t_contrast_change + rtLimit)); %attempting w/o lever_break column...
% corrIdxSum = nnz(corrIdx);
% IdxTable = table(corrIdx,EarlyBrkIdx,LtBrkIdx,EyeBrkIdx,PreTargLevBrk,PreCueLevBrk); % Table that has all the indexes, can make sure there are no overlaps or misses...
% %steps I used to determine that each row of the table only had one value,
% %and then found the index of rows that had more than one:
% % table2array(IdxTable)
% % sum(ans,2)
% %find(ans==2) << should have been 'find(ans>1)' but in any case, it worked.
% %Remember that 'ans' here are the answers for each iteration. I wasn't
% %storing anything as variables...
% 
% %6/5/23
% %DOING A CHECK ON THESE INDICES, I STILL END UP WITH AN OVERLAP SOMEWHERE
% %FOR 12 TRIALS (i.e. THERE ARE 12 MORE TRIALS THAN WHAT ARE LISTED FOR VALID)....NEED TO CONTINUE FIGURING OUT HOW TO FULLY ACCOUNT FOR
% %ALL TRIALS. NEED TO CHECK EACH ROW OF THE IdxTable TO SEE IF MORE THAN ONE
% %COLUMN CONTAINS 2 TRUES
% 
% 
% 
% 
% % ^^^ DIFFERENCE BETWEEN t_release AND t_lever_break ... trying to see what
% % these were...will come back to this.
% 
% 
% 
% 
% corrVal = Valid(corrIdx,:);
% PropCorrValid = height(corrVal)/(height(Valid)-(EyeBrksum+PreCueLevBrkSum)); % need to add in segmenting for times that aren't within the window but lever break still occurred, 
% 
% 
% %can use syntax like this:
% %'tLeverBreakVersion.t_contrast_change(12):tLeverBreakVersion.t_contrast_change(12)+150'
% 
% %ismember(Valid.break_eye,'false') & (Valid.t_lever_break < (Valid.t_contrast_change+150)) | ismember(Valid.break_eye,'false') & (Valid.t_lever_break > (Valid.t_contrast_change + rtLimit))
% 
% % Performance specific to each quadrant
% % example syntx: Valid.cue_x_loc == UL(1) & Valid.cue_y_loc == UL(2)
% 
% % There are a couple different types of trials which I am clueless as to
% % how they are happening. I think they are both still correct, but require
% % different logical operations in order to uncover...
% % A larger group seem to have the 'break_lever' column as set to true. In
% % these instances, t.lever_break == t.release. The release are also between
% % 150 - 800 ms after contrast change
% % There is also a second group though, whose 'break_lever' column is set to
% % false. In these instances 
% 
% 
% 
% 
% 
% 
% 
% % INVALID 
% 
% % Aggregate performance across all quadrants
% 
% % Performance specific to each quadrant
% 
% 
% % CATCH 
% 
% % Aggregate performance across all quadrants
% 
% % Performance specific to each quadrant
% 
%
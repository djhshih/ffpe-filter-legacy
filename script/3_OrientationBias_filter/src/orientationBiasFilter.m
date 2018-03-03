function [] = orientationBiasFilter(mafFilename, outputMAFFilename, outputDir, ...
    isCheckForMatFile, isGeneratingPlots, globalPbias, artifactThresholdRate, ...
    lod0Thresh, biasQP1,biasQP2,reference_allele, artifact_allele, bias_field)
% [] = orientationBiasFilter(mafFilename, outputMAFFilename, outputDir, isCheckForMatFile, isGeneratingPlots, globalPbias, artifactThresholdRate, lod0Thresh, biasQP1,biasQP2,reference_allele, artifact_allele, bias_field)
%
%   Entry point for running D-ToxoG generalized for any pair orientation bias.  
%       The generalization incarnation of the C>A/G>T artifact filter for any REF>ALT bias.
%       
%       NOTE: much of the code refers to the oxo-G artifact, but this
%       filter was rejiggered to deal with any orientation bias of the
%       oxo-G type: ie. ref>alt should be on F2R1 pairs, while the complement
%       SNVs appear on F1R2 pairs - governed by the binomal probability globalPbias
%
%   This script will perform the filtering and produce the following:
%   	MAF Files:  Filtered and Unfiltered
%           These will be named as <outputMAFFilename> and
%           <outputMAFFilename>.reject.maf, respectively.
%       Count Files:  Only contain the values of the total number of pass
%           and rejected mutations for the input maf file.
%           <outputMAFFilename>.pass_count.txt
%           <outputMAFFilename>.reject_count.txt
%           
%       Figures
%       Tables
%
%       The latter two can be used for further reporting downstream.
%
%   Note:  This script will preserve the order of mutations with the
%       following exception:
%
%       If multiple cases are included in the input maf file and the
%       mutations of cases are interleaved.  In this case, this script will
%       produce maf files with the mutations appearing clustered by case.
%
%       Simple example: 
%           Case input (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%
%           Case output (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
% 
% mafFilename -- input maf filename to be filtered.  This file will not be
%   modified.
%   Required Headers:
%     Chromosome -- Contig number without a prefix.  E.g. "3" or "X" (without quotes)
%     Start_position -- Position of the SNV
%     End_position -- Should be the same as Start_position for SNV
%     Reference_Allele -- The allele found in the reference genome at this position.  Single character representing the base, e.g. "G".
%     Tumor_Seq_Allele2 -- alternate allele.  Single character representing the base, e.g. "G".
%     Tumor_Sample_Barcode -- name of the tumor (or "case") sample.  This is used to generate file names and plots.
%     Matched_Norm_Sample_Barcode -- name of the normal (or "control") sample.  This is used to generate file names and plots.
%     ref_context -- Small window into the reference at the SNV.  The center position should be the same as Reference_Allele.  The total string should be of odd length and have a minimum length of 3. 
%       For example: Reference Allele is G, Chromosome is 1, Start_position and End_position are 120906037:  ref_context is CTTTTTTCGCGCAAAAATGCC  (string size is 21, in this case)
%     i_t_ALT_F1R2 -- the number of reads with pair orientation of F1R2 and with the alternate allele (Tumor_Seq_Allele2).
%     i_t_ALT_F2R1 -- the number of reads with pair orientation of F2R1 and with the alternate allele (Tumor_Seq_Allele2).
%     i_t_REF_F1R2 -- the number of reads with pair orientation of F1R2 and with the reference allele (Reference_Allele).
%     i_t_REF_F2R1 -- the number of reads with pair orientation of F2R1 and with the reference allele (Reference_Allele).
%         C>anything:  numerator is i_t_ALT_F2R1
%         A>anything:  numerator is i_t_ALT_F2R1
%         G>anything:  numerator is i_t_ALT_F1R2
%         T>anything:  numerator is i_t_ALT_F1R2
%     Variant_Type -- "SNP" (without quotes)
%
% outputMAFFilename -- filtered maf filename to be generated.  Without path
% information.
%
% outputDir (default: './') -- outputDir for filter results that are 
%   used for the report.
%   Note:  outputMAFFilename, figures, and tables will be located 
%       in outputDir.  The mat files will not be located in the outputDir.
% 
% globalPbias (default: .96) -- Pbias value to use in all estimations of
%   the given mafFile.  
%
% artifactThresholdRate (default .01) -- Expected proportion of mutations
%   that are artifacts.  e.g. .01 --> 1%  
%
% isCheckForMatFile (default: 0) -- looks for a mat file that contains the
%   mafTable already.  This is just to save on load time.  If value = 1, 
%   check for a mat file with name  mat/<mafFilename>.mat to load mafTable.
%   This option can be used to speed the running of this filter.
% 
% isGeneratingPlots (default: 0) -- true/false for whether this script
%   should generate the orchestra and lego png files.  These will be put in
%   the [outputDir]/figures/ .
%
% lod0Thresh (default: -1 --> No filtering) -- binomial mixture negative
%   log likelihood threshold.
%       The number being thresholded is the difference in the binomial
%       mixture model between the estimated weight of the oxoG component
%       and zero (the subtraction of the logs of each).  In other words, 
%       the log odds ratio of no oxoG contamination versus the oxoG
%       contamination (lod0).
%       If this difference is small, then this implies that there is not 
%       much OxoG contamination.  A large number implies lots of OxoG 
%       effect.  The difference can never be negative, since the oxoG 
%       estimate is the minimum.
%       Any cases where the difference is below this threshold are not
%           filtered.
%       The functionality around this parameter is not well tested.  
%       Use with caution.
%
%
% biasQP1 (default 1e8  - filter independent of biasQ) 
% biasQP2 (default 0   - filter independent of biasQ) 
%   Two fitting paramters to estimate the number of OxoG mutations from
%   biasQ value. biasQP1=1e8 and biasQP2=0.5 result in a filter that remains 
%   like the old filter below biasQ=35, becomes less
%   agressive above biasQ~35 and does no filtering above biasQ~45. 
%
%   biasQP2=0 forces a sharp cutoff for filtering only OxoG<biasQP1 samples.   
%
%
% reference_allele (default C)
% artifact_allele (default A)
%   reference and artifact allele expected to appear on F2R1 pairs while 
%   the complement appears on F1R2 pairs 
%
% bias_field (default 'picard_oxoQ')
%   prefix for all new fields Q,F,'_p','_q','_cut'
%
% See also load_table

[~, filename, ext] = fileparts(mafFilename);
if ~exist('mat', 'dir')
   mkdir('.', 'mat');
end
matFilename = ['mat/' filename ext '.mat'];

%% Parse arguments and determine which should use default values.

%if (nargin < 13 || ~exist('bias_field') || isempty(bias_field)), biasF_field = 'i_t_Foxog'; end
disp(['bias_field: ' bias_field ])

biasQ_field=[bias_field '_Q']
biasF_field=[bias_field '_F']
biasM_field=[bias_field '_mode']
biasp_field=[bias_field '_p_value']
biasq_field=[bias_field '_q_value']
biasC_field=[bias_field '_cut']


if (nargin < 12 || ~exist('reference_allele') || isempty(reference_allele)), reference_allele = 'C'; end
disp(['reference_allele: ' reference_allele ])

if (nargin < 11 || ~exist('artifact_allele') || isempty(artifact_allele)), reference_allele = 'A'; end
disp(['artifact_allele: ' artifact_allele ])

if (nargin < 10 || ~exist('biasQP2') || isempty(biasQP2)), biasQP2 = 0; end
if ischar(biasQP2),biasQP2 = str2num(biasQP2); end
disp(['biasQ P2: ' num2str(biasQP2) ])

if (nargin < 9 || ~exist('biasQP1') || isempty(biasQP1)), biasQP1 = 1e8; end
if ischar(biasQP1),biasQP1 = str2num(biasQP1); end
disp(['biasQ P1: ' num2str(biasQP1) ])

if nargin < 8 || ~exist('lod0Thresh') || isempty(lod0Thresh)
    lod0Thresh = -1;
end
if ischar(lod0Thresh)
    lod0Thresh = str2num(lod0Thresh);
end
disp(['lod0 Threshold: ' num2str(lod0Thresh) ])


if nargin < 7 || ~exist('artifactThresholdRate') || isempty(artifactThresholdRate)
    artifactThresholdRate = .01;
end
if ischar(artifactThresholdRate)
    artifactThresholdRate = str2num(artifactThresholdRate);
end

if nargin < 6 || ~exist('globalPbias') || isempty(globalPbias)
    globalPbias = .96;
end
if ischar(globalPbias)
    globalPbias = str2num(globalPbias);
end
    
if nargin < 5 || ~exist('isGeneratingPlots') || isempty(isGeneratingPlots)
    isGeneratingPlots = 0;
end

if nargin < 4 || ~exist('isCheckForMatFile') || isempty(isCheckForMatFile)
    isCheckForMatFile = 0;
end

if nargin < 3 || ~exist('outputDir') || isempty(outputDir)
    outputDir = './';
end

if ~exist(mafFilename, 'file') && (isCheckForMatFile && ~exist(matFilename, 'file'))
   error(['Given maf file does not exist: '  mafFilename])
end

%% Load the maf file
if ~isCheckForMatFile
    disp(['Loading ' mafFilename ' ...'])
    [mafTable] = loadMAFTable(mafFilename);
else
    if ~exist(matFilename, 'file')
        disp(['Loading ' mafFilename ' ...'])
        [mafTable] = loadMAFTable(mafFilename);
        disp(['Saving ' matFilename ' ... '])
        save(matFilename, 'mafTable');
    end
    disp(['Loading ' matFilename ' ...'])
    load(matFilename);
end

%% Retrieve all cases (one entry for each unique case, not one entry per
%   mutation)
disp(['Retrieving cases for ' mafFilename ' ...'])
[pairs] = retrieveUniqueCaseControlBarcodes(mafTable);
disp(['Retrieved ' num2str(length(pairs)) ' cases'])

% Target rate for artifacts that escape filtering
fdrThresh = artifactThresholdRate;  

% Measured p used in the B_A(ac,p) distibution in the binomial mixture
% model.
Pbias = globalPbias;

% Alt Allele Counts to use
acs = [3:50];

% Null distribution (real mutation) binomial parameter, p
p0 = .5;

%% Create the output directories
if ~strcmp(outputDir(end), '/') 
    outputDir = [outputDir '/'];
end

if ~exist(outputDir, 'dir')
   mkdir(outputDir); 
end

outputFigureDir =[ outputDir 'figures/'];
if ~exist(outputFigureDir, 'dir')
   mkdir('.', outputFigureDir); 
end

% Write the figure location to a file.
tmpFid = fopen([outputDir 'figureLoc.txt'], 'w');
fprintf(tmpFid, '%s', GetFullPath(outputFigureDir));
fclose(tmpFid);

outputTableDir = outputDir;
if ~exist(outputTableDir, 'dir')
   mkdir('.', outputTableDir); 
end

%% Create a file for storing the lego information
legoFidBefore=fopen([outputDir 'legoCasesBefore.txt'], 'w');
legoFid=fopen([outputDir 'legoCases.txt'], 'w');
if ~isGeneratingPlots
   disp('Lego raw data files will be empty when not generating plots.') 
end

%% Initialize and write a file for the parameters as specified
paramsFilename = [outputTableDir 'params.tsv'];
pfid = fopen( paramsFilename, 'w');
pHeaders = [{[bias_field '_Pbias']};{'artifactThreshold'};{'lod0Threshold'}];
fprintf(pfid, [strjoin('\t', pHeaders) '\n'])
fprintf(pfid,'%0.3f\t%0.2f\t%0.6f\n', Pbias, fdrThresh, lod0Thresh);
fclose(pfid)

%% Initialize a file to list the cases in the input maf
caseFid = fopen([outputDir 'allCases.txt'], 'w');

%% Initialize a table for summary data where each case is a row in a table.
tableFilename = [outputTableDir 'caseTableData.tsv'];
tfid = fopen( tableFilename, 'w');
headers = [{'case'}; {'control'}; {'N'}; {'N_A'}; {'N_nA'};{'N_filtered'};{'N_oxo'}; {'N_oxo_CI_low'}; {'N_oxo_CI_high'};{'N_mut'}; {'RR'}; {'RR_CI_low'}; {'RR_CI_hi'}; {'PA'}; {'PA_CI_low'}; {'PA_CI_hi'};{'Pbias'};{'Pbias_Est_CI_low'};{'Pbias_Est_CI_hi'};{'Proportion_N_filt'};{'lod0'};{'oxoQ'}];
headers=regexprep(headers, 'oxo',regexprep(bias_field,'^i_',''));

fprintf(tfid, [strjoin('\t', headers) '\n'])

totalRejectCount = 0;
totalPassCount = 0;

%% Initialize a map for all of the mutations
% This map will be used to add new columns to the maf files.
mutationMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

%% Initialize an array for storage of Passing Artifact & Rejected Real counts (and hi/low CIs) in a summary view
allPAs = zeros(length(pairs),3);
allRRs = zeros(length(pairs),3);
Ns = zeros(length(pairs),1);

%% For each unique case, run the filter
for i = 1:length(pairs)

    pair = pairs(i);    
    unameIndices = retrieveMutationsForGivenCase(pair.case, pair.control, mafTable);
    unameMafTable = trimStruct(mafTable, unameIndices);
    
    caseName = [pair.case ' v. ' pair.control];
    fprintf(caseFid, [caseName '\n']);
    
    disp(['Case ' num2str(i) ': ' caseName])
    
    % Create a list of the mutations segregated by artifact mode
    %  M_A are the artifact mode mutations
    %  M_nA are the non-artifact mode mutations
    
    [M, M_A, M_nA] = generateMutations_by_mode(mafTable, unameIndices, reference_allele, artifact_allele, biasQ_field, biasF_field, biasM_field);
        
    % fetch oxoQ value from first mutation in MAF for this sample
    biasQ=M.(biasQ_field)(1);    
    % estimate expected number of oxoG mutations from oxoQ fit 
    %eQ=oxoQP1*exp(-1*oxoQP2*oxoQ);
    % cut suppression factor based on oxoQ
    fQ=1./(1+exp(biasQP2*(biasQ-biasQP1)));
    % if oxoQP2 is not set, then don't correct fdr in cutLine
    if (biasQP2==0)
        fQ=1;
    end
    % Create a cut using False Discovery Rate 
    fdrStructs(i) = cutLineFDRbias(M, fdrThresh, Pbias, fQ, bias_field);
    
    %% Estimate num oxo artifacts
    [Noxo,Nmut,alf,alfci,NoxoCI, lod0] = oxogBinomialMixture(M_A.alt_read_count, round(M_A.(biasF_field) .* M_A.alt_read_count), acs, Pbias, p0);

    %% If we are 95% certain that this case has less than 1% artifacts, do not filter.  Do this by setting 
    isFilterThisCase = 1;
    if lod0 < lod0Thresh
        disp('The logarithmic odds ratio is not large enough to warrant filtering this case.  Doing no filtering.')
        isFilterThisCase = 0;
    end
    
    if (biasQP2==0)&(biasQ>biasQP1)        
        disp(sprintf('biasQ exceeds threshold: %f>%f   Doing no filtering.',biasQ>biasQP1))
        isFilterThisCase = 0;
    end
    
    if NoxoCI(2) < (fdrThresh*length(M.alt_read_count))
       disp('We are 95% certain that this case has less than 1% artifacts.  Doing no filtering.')
       isFilterThisCase = 0;
    end
    if ~isFilterThisCase
       fdrStructs(i).(biasC_field) = zeros(size(fdrStructs(i).(biasC_field)));
    end
    
    %% Update pass and reject counts after finalizing fdrStructs(i).cut
    totalPassCount = totalPassCount + sum(~fdrStructs(i).(biasC_field));
    totalRejectCount = totalRejectCount + sum(fdrStructs(i).(biasC_field));

    
    %% Estimate the a cutline for FoxoG x alt count
    % acVal ./ acs is FoxoG.  The line is (FoxoG(i), acs(i))
    [acVal] = estimateCutLineFromFdrOutputBias(fdrStructs(i), acs, Pbias, biasp_field, biasC_field);
    cutLineBias = (acVal' ./ acs);

    %% Plotting figures
    if isGeneratingPlots
        
        %% Plot the artifact mode mutations
        disp('Plotting artifact mode (with patch)...')
        plotArtifactModeBias(outputFigureDir, caseName, M_A, fdrStructs(i).(biasC_field), cutLineBias, acs,biasQ_field, biasF_field)

        %% Plot the non-artifact mode mutations
        disp('Plotting non-artifact mode ...')
        plotArtifactModeBias(outputFigureDir, caseName, M_nA, [],[],acs,biasQ_field, biasF_field)

        %% Generate Lego plots (before and after)
        disp('Making lego plots...')
        Pbefore.label = caseName;
        figure('visible','off');
        legoStructsBefore = plotMutationSpectrumCategLego1(unameMafTable,'unit',[outputFigureDir caseName '_lego_before'], Pbefore);        
        close(gcf);
        saveLegoStructToFile(caseName, legoStructsBefore, legoFidBefore, i==1);
        
        Pafter.label = [caseName ' Filtered'];
        figure('visible','off');
        legoStructsAfter = plotMutationSpectrumCategLego1(trimStruct(unameMafTable, ~fdrStructs(i).(biasC_field)),'unit',[outputFigureDir caseName '_lego_after'], Pafter);
        close(gcf);
        saveLegoStructToFile(caseName, legoStructsAfter, legoFid, i==1);
    else
        disp('Figures are not being generated.')
    end
    
    %% Estimate false positive (i.e. rejected real)
    [FP,FPci,FN,FNci] = estimateFalseRatesForOxoGfdrBias(M.alt_read_count, M.(biasM_field), fdrStructs(i), Noxo, NoxoCI, Pbias,biasp_field,biasC_field);
    allPAs(i,:) = [FN FNci(1) FNci(2)]; 
    allRRs(i,:) = [FP FPci(1) FPci(2)]; 
    Ns(i) = length(M.alt_read_count);
    
    %% Estimate Pbias error.  
    % Note that this estimate is not actually used.  This estimate is
    %  recorded for later analysis.
    [PbiasEst,PbiasEstci] = PoxoGestimate(M.alt_read_count, round(M.(biasF_field) .* M.alt_read_count),fdrStructs(i).(biasM_field), [], Noxo, NoxoCI);
    
    %% Write stdout
    tbl_entry = [strjoin('\t', [{pair.case}; {pair.control};{num2str(length(M.(biasF_field)))}; {num2str(length(M_A.(biasF_field)))}; ...
        {num2str(length(M_nA.(biasF_field)))}; {num2str( sum(fdrStructs(i).(biasC_field)))}; {num2str(Noxo)}; {num2str(NoxoCI(1))}; {num2str(NoxoCI(2))}; {num2str(Nmut)}; ...
        {num2str(FP)}; {num2str(FPci(1))}; {num2str(FPci(2))};{num2str(FN)}; {num2str(FNci(1))}; {num2str(FNci(2))}; ...
        {num2str(PbiasEst)}; {num2str(PbiasEstci(1))}; {num2str(PbiasEstci(2))} ;{num2str( sum(fdrStructs(i).(biasC_field))/length(M.(biasF_field)) )};{num2str(lod0)};...
        {num2str(biasQ)}...
        ]) '\n'];
    fprintf(1,  tbl_entry)
    
    %% Write table entry to table file
    fprintf(tfid, tbl_entry )
    
    %% Generate temporary table for later use in appending (column-wise) to output maf files
    caseMutMap = addMutationEntriesBias(pair, unameMafTable, M, fdrStructs(i), acVal, acs, Pbias,biasp_field,biasq_field,biasC_field,biasM_field,biasF_field);
    mutationMap = [mutationMap; caseMutMap];
end
fclose(legoFidBefore);
fclose(legoFid);
fclose(tfid);
fclose(caseFid);

%% Write pass and reject counts to a file.
writeFileWithSingleNumber([outputDir  outputMAFFilename '.pass_count.txt'], totalPassCount);
writeFileWithSingleNumber([outputDir  outputMAFFilename '.reject_count.txt'], totalRejectCount);

%% Write out new maf files (one for pass and one for reject)
%   Additionally, writeFilterMAFFiles appends new columns to the mafFile,
%       so overwrite the mafTable variable to include these.
disp('Writing MAF files and populating table structure with new columns...')
mafTable = writeFilterMAFFilesBias(mafTable, [outputDir outputMAFFilename], mutationMap,biasp_field,biasq_field,biasC_field,biasM_field,biasF_field);

stub = regexprep(bias_field,'^i_','');


%% Plot summary data across all cases
if (isGeneratingPlots) && (length(pairs) > 0)
    disp(' Generating summary plots for all cases...')
    [dummy, basename, ext] = fileparts(outputMAFFilename);
    plotFalseRates([outputFigureDir outputMAFFilename '.RR'], [basename ' Rejected Reals'], 'Rejected Real Count', Ns, allRRs);
    plotFalseRates([outputFigureDir outputMAFFilename '.PA'], [basename ' Passing Artifacts'], 'Passing Artifact Count', Ns, allPAs);
    disp('False rate plots completed.')
    
    % Plot summary orchestra plot of all artifact mode mutations
    allFdrStructs = fdrStructs(1);
    for k = 2:length(fdrStructs)
        allFdrStructs = mergeStruct(allFdrStructs, fdrStructs(k));
    end
    [M_summary, M_A_summary, M_nA_summary] = generateMutationsBias(mafTable, [1:length(mafTable.Start_position)], reference_allele, artifact_allele,biasQ_field,biasM_field,biasF_field);
    plotArtifactModeBias(outputFigureDir, [outputMAFFilename ' All Artifact Mode'], M_A_summary, allFdrStructs.(biasC_field) , [], acs, biasQ_field, biasF_field)
    plotArtifactModeBias(outputFigureDir, [outputMAFFilename ' All Non-Artifact Mode'], M_nA_summary, [],[],acs,biasQ_field, biasF_field)
    disp('Summary Orchestra plots completed.')
    
    Pbefore.label = basename;
    figure('visible','off');
    [legoStructBeforeSet] = plotMutationSpectrumCategLego1(mafTable,'unit',[outputFigureDir outputMAFFilename '_lego_before'], Pbefore);
    close(gcf);
    Pafter.label = [basename ' Filtered'];
    figure('visible','off');
    [legoStructAfterSet] = plotMutationSpectrumCategLego1(trimStruct(mafTable, find(~mafTable.(biasC_field))),'unit',[outputFigureDir outputMAFFilename '_lego_after'], Pafter)
    close(gcf);
    disp('Summary Lego plots completed.')
    
    summaryTable = load_table(GetFullPath(tableFilename), char(9));
    figure('visible','off');
    hist(summaryTable.Proportion_N_filt, [.005:.01:.995])
    xlabel('N_Filtered/N', 'FontSize', 14, 'Interpreter', 'None')
    ylabel('Count', 'FontSize', 14, 'Interpreter', 'None')
    title(['Proportion of Filtered Mutations'], 'FontSize', 16, 'Interpreter', 'None')
    print('-dpng', [outputFigureDir outputMAFFilename '_filteredProportion_hist.png']);
    close(gcf)
    disp('Proportion of Filtered Mutations plot completed.')
    
    % filter summary
    kpass0=find(mafTable.(biasq_field)<fdrThresh);
    kpass=find(~mafTable.(biasC_field));
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode));
    t=[];
    t.sample=tmp(:,1);
    t.ntot=cell2mat(tmp(:,2));
    [i m]=ismember(t.sample,mafTable.Tumor_Sample_Barcode);
    t.(biasQ_field)=NaN*t.ntot;
    t.(biasQ_field)(i,1)=mafTable.(biasQ_field)(m(m>0));
    i=t.(biasQ_field)>50; % spurious oxoQ when error rate is negative
    t.(biasQ_field)(i)=rand(sum(i),1)*0.5+59.5;
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode(kpass0)));
    [i m]=ismember(t.sample,tmp(:,1));
    t.npass0=0*t.ntot;
    t.npass0(i)=cell2mat(tmp(m(m>0),2));
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode(kpass)));
    [i m]=ismember(t.sample,tmp(:,1));
    t.npass=0*t.ntot;
    t.npass(i)=cell2mat(tmp(m(m>0),2));
    figure('visible','off')
    subplot('position',[0.1 0.53 0.85 0.39])
    h=semilogy(t.(biasQ_field),t.ntot*0,'go',t.(biasQ_field),t.ntot-t.npass0+1,'r+',t.(biasQ_field),t.ntot-t.npass+1,'bx','markersize',8,'linewidth',2);
    set(h(1),'color',0.5*[0.1 1 0.1])
    set(gca,'xticklabel',[])
    ylabel('Mutations cut +1', 'FontSize', 14, 'Interpreter', 'None')
    stub=regexprep(bias_field,'^i_','');
    legend(['No ' stub ' filter'],['filter on ' stub ' q-value only'],['final filter using ' stub],'location','NE')
    title(['Filter dependence on ' biasQ_field], 'FontSize', 16, 'Interpreter', 'None')
    xlim([0 60])
    grid on
    subplot('position',[0.1 0.12 0.85 0.39])
    h=semilogy(t.(biasQ_field),t.ntot+1,'go',t.(biasQ_field),t.npass0+1,'r+',t.(biasQ_field),t.npass+1,'bx','markersize',8,'linewidth',2);
    set(h(1),'color',0.5*[0.1 1 0.1])
    xlabel(biasQ_field, 'FontSize', 14, 'Interpreter', 'None')    
    ylabel('Mutations passed +1', 'FontSize', 14, 'Interpreter', 'None')
%     k=find(xT>40); 
%     if ~isempty(k)
%         xTL{k(end)}='~ 80';
%         set(gca,'xticklabel',xTL);
%     end
    xlim([0 60])    
    xL=get(gca,'xlim');    xT=get(gca,'xtick');    xTL=cellstr(get(gca,'xticklabel')); 

    ylim([min([10; t.npass/2]) max([2e3; t.npass*2])])
    yL=get(gca,'ylim');
    line(40*[1 1],yL,'color',0.5*[.5 0.5 1],'linestyle',':')
    line(xL,xL*0+median(t.npass),'color',0.5*[.5 0.5 1],'linestyle','--')
    text(0.99,0.95,sprintf(' median mutations passed = %d ',round(median(t.npass))),...
        'units','normalized', 'FontSize', 14, 'HorizontalAlignment','right',... 
        'BackgroundColor','w','edgecolor','k')
  
    grid on
    print('-dpng', [outputFigureDir outputMAFFilename '.' biasQ_field '.png']);
    close(gcf)
    disp([biasQ_field ' plot completed.'])
    
    % allele frequency
    plotAllelicFractionSummaries(outputFigureDir, basename, [outputMAFFilename '_AF'], mafTable,biasC_field)
    disp('Summary AF plots completed.')
    

    
end

%% Clean up
close all    % close all figures so xvnc can terminate
disp('Done.')

    

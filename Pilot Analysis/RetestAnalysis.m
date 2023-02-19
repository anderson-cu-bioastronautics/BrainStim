function [RetestData] = RetestAnalysis(subID,TestType,SRType,optimal)

%% Generate the file names to analyze
currentFolder = pwd;
typeFolder = 'RA1Data';
subFold = num2str(subID);

switch TestType
    case 1 % define visual folder information
        mu_guess = 0.2; sigma_guess = 0.05; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Visual_0mA_'];
            parseFold2=['dB.xls'];
            parseFold3=['dB_2.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Visual_'];
            parseFold2=['mA_0dB.xls'];
            parseFold3=['mA_0dB_2.xls'];
        else
        end
    case 2 % define auditory folder information
        mu_guess = 8; sigma_guess = 3; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Auditory_0mA_'];
            parseFold2=['dB.xls'];
            parseFold3=['dB_2.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Auditory_'];
            parseFold2=['mA_0dB.xls'];
            parseFold3=['mA_0dB_2.xls'];
        else
        end
    case 3 % define tactile folder information
        mu_guess = 0.3; sigma_guess = 0.1; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Tactile_0mA_'];
            parseFold2=['dB.xls'];
            parseFold3=['dB_2.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Tactile_'];
            parseFold2=['mA_0dB.xls'];
            parseFold3=['mA_0dB_2.xls'];
        else
        end
    case 4
        mu_guess = 0; sigma_guess = 0; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Vestibular_0mA_'];
            parseFold2=['dB.xls'];
            parseFold3=['dB_2.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Vestibular_'];
            parseFold2=['mA_0dB.xls'];
            parseFold3=['mA_0dB_2.xls'];
        else
        end
end

% Original sham
sham1fid = [parseFold1,'0',parseFold2];

% Retest sham
sham2fid = [parseFold1,'0',parseFold3];

% Original optimal
if SRType == 2 && optimal == 1
    opt1fid = [parseFold1,num2str(1),parseFold2];
elseif sign(optimal) == -1
    opt1fid = [parseFold1,'neg',num2str(abs(optimal)),parseFold2];
elseif floor(optimal)~= optimal
    opt1fid = [parseFold1,num2str(floor(optimal)),'p',num2str(10*rem(optimal,1)),parseFold2];
else
    opt1fid = [parseFold1,num2str(abs(optimal)),parseFold2];
end

% Retest optimal
if SRType == 2 && optimal == 1
    opt2fid = [parseFold1,num2str(1),parseFold3];
elseif sign(optimal) == -1
    opt2fid = [parseFold1,'neg',num2str(abs(optimal)),parseFold3];
elseif floor(optimal)~= optimal
    opt2fid = [parseFold1,num2str(floor(optimal)),'p',num2str(10*rem(optimal,1)),parseFold3];
else
    opt2fid = [parseFold1,num2str(abs(optimal)),parseFold3];
end

%% Read in the data

% Original sham
sham1file = fullfile(currentFolder,typeFolder,subFold,sham1fid);
sham1master = importdata(sham1file);
sham1stim = sham1master.data.Sheet1(:,2);
if TestType ~= 2
    sham1stim = abs(sham1stim);
end
sham1cor = sham1master.data.Sheet1(:,4);

% Retest sham
sham2file = fullfile(currentFolder,typeFolder,subFold,sham2fid);
sham2master = importdata(sham2file);
sham2stim = sham2master.data.Sheet1(:,2);
if TestType ~= 2
    sham2stim = abs(sham2stim);
end
sham2cor = sham2master.data.Sheet1(:,4);

% Pooled sham
shamStimPooled = [sham1stim;sham2stim];
shamCorPooled = [sham1cor;sham2cor];

% Original optimal
opt1file = fullfile(currentFolder,typeFolder,subFold,opt1fid);
opt1master = importdata(opt1file);
opt1stim = opt1master.data.Sheet1(:,2);
if TestType ~= 2
    opt1stim = abs(opt1stim);
end
opt1cor = opt1master.data.Sheet1(:,4);

% Retest optimal
opt2file = fullfile(currentFolder,typeFolder,subFold,opt2fid);
opt2master = importdata(opt2file);
opt2stim = opt2master.data.Sheet1(:,2);
if TestType ~= 2
    opt2stim = abs(opt2stim);
end
opt2cor = opt2master.data.Sheet1(:,4);

% Pooled optimal
optStimPooled = [opt1stim;opt2stim];
optCorPooled = [opt1cor;opt2cor];

%% Run the T-tests

% Guesses
guesses = [mu_guess sigma_guess 0];
lambda_max = 0.06;

% Original sham and retest sham
t_test_sham1_sham2 = t_test_satterthwaite(sham1stim,sham2stim,sham1cor,sham2cor,guesses,lambda_max);

% Original optimal and retest optimal
t_test_opt1_opt2 = t_test_satterthwaite(opt1stim,opt2stim,opt1cor,opt2cor,guesses,lambda_max);

% Original sham and original optimal
t_test_sham1_opt1 = t_test_satterthwaite(sham1stim,opt1stim,sham1cor,opt1cor,guesses,lambda_max);

% Original sham and retest optimal
t_test_sham1_opt2 = t_test_satterthwaite(sham1stim,opt2stim,sham1cor,opt2cor,guesses,lambda_max);

% Retest sham original optimal
t_test_sham2_opt1 = t_test_satterthwaite(sham2stim,opt1stim,sham2cor,opt1cor,guesses,lambda_max);

% Retest sham retest optimal
t_test_sham2_opt2 = t_test_satterthwaite(sham2stim,opt2stim,sham2cor,opt2cor,guesses,lambda_max);

% Retest sham retest optimal
t_test_pooled = t_test_satterthwaite(shamStimPooled,optStimPooled,shamCorPooled,optCorPooled,guesses,lambda_max);

RetestData = struct();
RetestData.sham1sham2 = t_test_sham1_sham2;
RetestData.opt1opt2 = t_test_opt1_opt2;
RetestData.sham1opt1 = t_test_sham1_opt1;
RetestData.sham2opt1 = t_test_sham2_opt1;
RetestData.sham1opt2 = t_test_sham1_opt2;
RetestData.sham2opt2 = t_test_sham2_opt2;
RetestData.pooledSham_pooledOpt = t_test_pooled;
end
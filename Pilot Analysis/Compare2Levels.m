% Function to analyze and observe two threshold conditions

function []=Compare2Levels(subID,TestType,SRType,levels)

%% Observe singular file data and make error bars

currentFolder=pwd;
typeFolder='RA1Data';
subFold=num2str(subID);

%% Find the noise level that was tested twice
switch TestType
    case 1 % define visual folder information
        mu_guess = 0.2; sigma_guess = 0.05; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Visual_0mA_'];
            parseFold2=['*dB_2.xls'];
            label2=['dB.xls'];
            label2a=['dB_2.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Visual_'];
            parseFold2=['*mA_0dB_2.xls'];
            label2=['mA_0dB.xls'];
            label2a=['mA_0dB_2.xls'];            
        else
        end
    case 2 % define auditory folder information
        mu_guess = 8; sigma_guess = 3; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Auditory_0mA_'];
            parseFold2=['*dB_2.xls'];
            label2=['dB.xls'];
            label2a=['dB_2.xls'];            
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Auditory_'];
            parseFold2=['*mA_0dB_2.xls'];
            label2=['mA_0dB.xls'];
            label2a=['mA_0dB_2.xls'];            
        else
        end
    case 3 % define tactile folder information
        mu_guess = 0.3; sigma_guess = 0.1; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Tactile_0mA_'];
            parseFold2=['*dB_2.xls'];
            label2=['dB.xls'];
            label2a=['dB_2.xls'];            
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Tactile_'];
            parseFold2=['*mA_0dB_2.xls'];
            label2=['mA_0dB.xls'];
            label2a=['mA_0dB_2.xls'];            
        else
        end
    case 4
        mu_guess = 0; sigma_guess = 0; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Vestibular_0mA_'];
            parseFold2=['*dB_2.xls'];
            label2=['dB.xls'];
            label2a=['dB_2.xls'];            
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Vestibular_'];
            parseFold2=['*mA_0dB_2.xls'];
            label2=['mA_0dB.xls'];
            label2a=['mA_0dB_2.xls'];
        else
        end
end

% Identify the folders appended with _2
checkFile = [parseFold1 parseFold2];
twoLevel = dir(fullfile(currentFolder,typeFolder,subFold,checkFile));
TestFilename = getfield(twoLevel,{2},'name');
strTestLevel = TestFilename; % Pull out noise level
strTestLevel(1:length(parseFold1))=[];
strTestLevel(end-(length(parseFold2)-2):end)=[];

if length(strTestLevel) == 3 % means it's a decimal
   decLevel = ['0.' strTestLevel(end)];
   secondLevel = str2num(decLevel); % get actual level
   secondLevel = round(secondLevel,1);
else
   secondLevel = str2num(strTestLevel);
end
bothLevel=[0 secondLevel];
bothLevelstr={[num2str(bothLevel(1)) 'a']; [num2str(bothLevel(1)) 'b']; [num2str(bothLevel(2)) 'a']; [num2str(bothLevel(2)) 'b']};
combLevelstr={[num2str(bothLevel(1)) 'comb'] [num2str(bothLevel(2)) 'comb']};
%% Find thresholds and organize individual data, run statistics against each other

for i=1:2
    level=bothLevel(i);
 if sign(level) == -1 % folder definitions for negative SR levels
        fileID{i}=[parseFold1 'neg' num2str(abs(level)) label2];
        fileIDa{i}=[parseFold1 'neg' num2str(abs(level)) label2a];
    elseif floor(level) ~= level % definitions for non-integer SR levels
        fileID{i}=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) label2];
        fileIDa{i}=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) label2a];
    else % definition for all whole number SR levels
        fileID{i}=[parseFold1 num2str(abs(level)) label2];
        fileIDa{i}=[parseFold1 num2str(abs(level)) label2a];        
 end
end

fid1=fullfile(currentFolder,typeFolder,subFold,fileID{1});
fid1a=fullfile(currentFolder,typeFolder,subFold,fileIDa{1});
fid2=fullfile(currentFolder,typeFolder,subFold,fileID{2});
fid2a=fullfile(currentFolder,typeFolder,subFold,fileIDa{2});
       

[X1,cor1,mu_bi1,sigma_bi1,mu_hat1,sigma_hat1,x1]=SingleFileAnalysis(fid1,mu_guess,sigma_guess,subID,TestType);
[X2,cor2,mu_bi2,sigma_bi2,mu_hat2,sigma_hat2,x2]=SingleFileAnalysis(fid2,mu_guess,sigma_guess,subID,TestType);

[X1a,cor1a,mu_bi1a,sigma_bi1a,mu_hat1a,sigma_hat1a,x1a,]=SingleFileAnalysis(fid1a,mu_guess,sigma_guess,subID,TestType);
[X2a,cor2a,mu_bi2a,sigma_bi2a,mu_hat2a,sigma_hat2a,x2a,]=SingleFileAnalysis(fid2a,mu_guess,sigma_guess,subID,TestType);

% arrange our data
X=[X1 X1a X2 X2a];
cor=[cor1 cor1a cor2 cor2a];
mu_bi=[mu_bi1 mu_bi1a mu_bi2 mu_bi2a];
sigma_bi=[sigma_bi1 sigma_bi1a sigma_bi2 sigma_bi2a];
mu_hat=[mu_hat1 mu_hat1a mu_hat2 mu_hat2a];
sigma_hat=[sigma_hat1 sigma_hat1a sigma_hat2 sigma_hat2a];

% cycle through jack knife for binary and only binary
for i=1:4
musi_bi = jackknife(@just794_bi,X(:,i),cor(:,i),mu_guess,sigma_guess);
mu_se_bi(i) = sqrt( (N-1)/N*sum( (mu_bi(i) - musi_bi(:,1)).^2) );
mu_jackData{i} = musi_bi(:,1);
mu_indDOF(i)=sum((mu_bi(i) - musi_bi(:,1)).^2)^2/sum((mu_bi(i) - musi_bi(:,1)).^4);
sigma_se_bi(i) = sqrt( (N-1)/N*sum( (sigma_bi(i) - musi_bi(:,2)).^2 ) );
end

% do some statistical comparisons between the 2 individuals and against
% each other for the individuals
count=2;
tab_cond{1,1}='File 1'; tab_cond{1,2}='File 2'; tab_cond{1,3}='t value'; tab_cond{1,4}='p value'; tab_cond{1,5}='DOF';
       
for i=1:4
    if i~=4
   for j=i+1:4
       mu_compDOF(count)=(mu_se_bi(i)^2+mu_se_bi(j)^2)^2/(1/mu_indDOF(i)*mu_se_bi(i)^4+1/mu_indDOF(j)*mu_se_bi(j)^4); % real DOF
       t_ind(count)=(mu_bi(i)-mu_bi(j))/sqrt(mu_se_bi(i)^2+mu_se_bi(j)^2);
       p_ind(count)=1-tcdf(t_ind(count),mu_compDOF(count));
       tab_cond{count,1}=bothLevelstr{i}; tab_cond{count,2}=bothLevelstr{j};
       tab_cond{count,3}=t_ind(count); tab_cond{count,4}=p_ind(count); tab_cond{count,5}=mu_compDOF(count);
       count=count+1;
   end
    end
end


%% Observe combined file data and make error bars

% Look at sham condition
[Xall1,corall1,mu_biall1,sigma_biall1,mu_hatall1,sigma_hatall1]=DoubleFileAnalysis(fid1,fid1a,mu_guess,sigma_guess,subID,TestType);

[musi_biall1] = jackknife(@just794_bi,Xall1,corall1,mu_guess,sigma_guess);
mu_se_biall1 = sqrt( (N-1)/N*sum( (mu_biall1 - musi_biall1(:,1)).^2) );
sigma_se_biall1 = sqrt( (N-1)/N*sum( (sigma_biall1 - musi_biall1(:,2)).^2 ) );
mu_DOFcom1=sum((mu_biall1 - musi_biall1(:,1)).^2)^2/sum((mu_biall1 - musi_biall1(:,1)).^4);

% Look at best condition
[Xall2,corall2,mu_biall2,sigma_biall2,mu_hatall2,sigma_hatall2]=DoubleFileAnalysis(fid2,fid2a,mu_guess,sigma_guess,subID,TestType);

[musi_biall2] = jackknife(@just794_bi,Xall2,corall2,mu_guess,sigma_guess);
mu_se_biall2 = sqrt( (N-1)/N*sum( (mu_biall2 - musi_biall2(:,1)).^2) );
sigma_se_biall2 = sqrt( (N-1)/N*sum( (sigma_biall2 - musi_biall2(:,2)).^2 ) );
mu_DOFcom2=sum((mu_biall2 - musi_biall2(:,1)).^2)^2/sum((mu_biall2 - musi_biall2(:,1)).^4);

% analyze and document combined file
mu_compDOF(count)=(mu_se_biall1^2+mu_se_biall2^2)^2/(1/mu_DOFcom1*mu_se_biall1^4+1/mu_DOFcom2*mu_se_biall2^4); % real DOF
t_ind(count)=(mu_biall1-mu_biall2)/sqrt(mu_se_biall1^2+mu_se_biall2^2);
p_ind(count)=1-tcdf(t_ind(count),mu_compDOF(count));
tab_cond{count,1}=combLevelstr{1}; tab_cond{count,2}=combLevelstr{2};
tab_cond{count,3}=t_ind(count); tab_cond{count,4}=p_ind(count); tab_cond{count,5}=mu_compDOF(count);

tab_cond



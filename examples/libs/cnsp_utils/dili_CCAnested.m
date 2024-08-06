function [bestDprimeSVM,bestRR,score] = diliCCAnested(stim,eeg,nCompStimAll,nCompEEGAll,shifts,prefixTmpFile)
    if ~exist('prefixTmpFile')
        prefixTmpFile = '';
    end
    % Precalculating covariance matrices
%     disp('Pre-computing covariance matrices')
%     for iEEG = 1:length(nCompEEGAll)
%         nCompEEG = nCompEEGAll(iEEG);
%         yy = eegDimReduction(eeg',nCompEEG);
%         save(['tmpCCAData',prefixTmpFile,'_yy_',num2str(iEEG),'.mat'],'yy')
%         
%         for iStim = 1:length(nCompStimAll)
%             disp(sprintf('\b.'))
%             nCompStim = nCompStimAll(iStim);
%             if iEEG == 1
%                 xx = stimDimReduction(stim,nCompStim);
%                 save(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx')
%             else
%                 load(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx')
%             end
%             
%             nTrials=length(eeg);
%             n=size(xx{1},2)+size(yy{1},2);
%             C=zeros(n,n,length(shifts),nTrials);
%             for iTrial=1:nTrials
%                 C(:,:,:,iTrial)=nt_cov_lags(xx{iTrial}, yy{iTrial}, shifts);
%             end
% 
%             save(['tmpCCAData',prefixTmpFile,'_C_',num2str(iStim),'_',num2str(iEEG),'.mat'],'C')
%         end
%     end
    
    % Tuning
    % Outer cross-validation loop
    if length(nCompEEGAll) == 1 && length(nCompStimAll) == 1
        autoTuning = 0;
    else
        autoTuning = 1;
    end
    
    clear nCompStimIdxBest nCompEEGIdxBest
%     for tr = 1:length(eeg)
    tr = 1;
        disp(sprintf('\b.'))

        trainTr = setdiff(1:length(eeg),tr);
        testTr = tr;

        % This is not a full search. It's searching only around the first
        % local optimum
        clear bestShiftAll score AA BB
        idxToExplore = 1:2:length(nCompStimAll); % Initially, explore only a few idx, fine tuning from the next step
        stimCovExist = zeros(length(nCompStimAll),1);
        for iEEG = 1:length(nCompEEGAll)
            nCompEEG = nCompEEGAll(iEEG);
%             load(['tmpCCAData',prefixTmpFile,'_yy_',num2str(iEEG),'.mat'],'yy')
            disp(sprintf('\b,'))
            yy = eegDimReduction(eeg',nCompEEG);
            if autoTuning, save(['tmpCCAData',prefixTmpFile,'_yy_',num2str(iEEG),'.mat'],'yy'); end
            
            for iStim = idxToExplore %1:length(nCompStimAll)
                disp(sprintf('\b;'))
                nCompStim = nCompStimAll(iStim);
%                 load(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx')
                if ~stimCovExist(iStim)
                    xx = stimDimReduction(stim,nCompStim);
                    stimCovExist(iStim) = 1;
                    if autoTuning, save(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx'); end
                else
                    load(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx')
                end

                % Getting cov matrix
%                 load(['tmpCCAData',prefixTmpFile,'_C_',num2str(iStim),'_',num2str(iEEG),'.mat'],'C')
                nTrials=length(eeg);
                n=size(xx{1},2)+size(yy{1},2);
                C=zeros(n,n,length(shifts),nTrials);
                for iTrial=1:nTrials
                    C(:,:,:,iTrial)=nt_cov_lags(xx{iTrial}, yy{iTrial}, shifts);
                end
                if autoTuning, save(['tmpCCAData',prefixTmpFile,'_C_',num2str(iStim),'_',num2str(iEEG),'.mat'],'C'); end

                % Determining best shift
                [AA{iStim,iEEG},BB{iStim,iEEG},RR] = ...
                    nt_cca_crossvalidate_Gio(xx(trainTr),yy(trainTr),-shifts,0,C(:,:,:,trainTr));
                meanRR = mean(RR,3); % mean across trials
    %             figure;plot(shifts/128*1000,sum(meanRR,1))
                [bestR,bestShiftAll(iStim,iEEG)] = max(sum(meanRR(1:5,:),1));
%                 [bestR,bestShiftAll(iStim,iEEG)] = max(sum(meanRR(1:10,:),1));
                bestShift = bestShiftAll(iStim,iEEG);
                disp(['bestShift = ',num2str(shifts(bestShift)/64*1000),' ms'])
                % dprime for training set
%                 [~,~,~,dprime,dprime2,dprimeSVM]=nt_cca_crossvalidate_Gio(xx(trainTr),yy(trainTr),-shifts(bestShift),1,C(:,:,:,trainTr));

%                 score(iStim,iEEG) = dprimeSVM(4)
                score(iStim,iEEG) = bestR

                stimRank = find(diff(sum(xx{1},1))==0,1);
                disp(['stimRank ',num2str(stimRank),'; nCompStim ',num2str(nCompStim)])
                if stimRank < nCompStim % we have reached the rank, no point in increasing nStim
                    break
                end
            end
            [~,maxEEGIdx] = max(score(:,iEEG));
            idxToExplore = [maxEEGIdx-1,maxEEGIdx,maxEEGIdx+1];
            idxToExplore = idxToExplore(idxToExplore>0 & idxToExplore<=length(nCompStimAll));
        end

    
    iStimPrev = 0;
    iEEGPrev = 0;
    for tr = 1:length(eeg) % keeping the same parameters (they won't change unless the detail is sufficient) - to change for final result
        trainTr = setdiff(1:length(eeg),tr);
        testTr = tr; % move the for loop up ahead in the final version / here we assume that the params don't change
        
        % Selecting best parameters
        maxValue = max(score(:));
        [nCompStimIdxBest(tr),nCompEEGIdxBest(tr)] = find(score == maxValue,1);
        iStim = nCompStimIdxBest(tr);
        iEEG = nCompEEGIdxBest(tr);
        bShift = -shifts(bestShiftAll(iStim,iEEG));
        if autoTuning
            if iStimPrev ~= iStim, load(['tmpCCAData',prefixTmpFile,'_xx_',num2str(iStim),'.mat'],'xx'); end
            if iEEGPrev ~= iEEG, load(['tmpCCAData',prefixTmpFile,'_yy_',num2str(iEEG),'.mat'],'yy'); end
            if iEEGPrev ~= iEEG || iStimPrev ~= iStim, load(['tmpCCAData',prefixTmpFile,'_C_',num2str(iStim),'_',num2str(iEEG),'.mat'],'C'); end
        end
        [bestRR{tr},bestDprime{tr},bestDprime2{tr},bestDprimeSVM{tr}] = ...
            nt_cca_test(xx{testTr},yy{testTr},xx(trainTr),yy(trainTr),bShift,1,C(:,:,:,trainTr));
%         bestDprimeSVM{tr}(4)
        
%         [bestRR{tr},bestDprime{tr},bestDprime2{tr},bestDprimeSVM{tr}] = ...
%             nt_cca_test(xx{testTr},yy{testTr},xx(trainTr),yy(trainTr),AA{iStim},BB{iEEG},bShift,1);
% %         [~,~,~,dprime,dprime2,dprimeSVM] = ...
%             nt_cca_crossvalidate(xx(trainTr),yy(trainTr),-shifts(bestShift),1);
        iStimPrev = iStim;
        iEEGPrev = iEEG;
    end
end

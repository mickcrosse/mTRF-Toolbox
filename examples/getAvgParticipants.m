% Cognition and Natural Sensory Processing (CNSP) Workshop
% Example 1 - Forward TRF
%
% This example script loads and preprocesses a publicly available dataset
% (you can use any of the dataset in the CNSP resources). Then, the script
% runs a typical forward TRF analysis.
%
% Note:
% This code was written with the assumption that all subjects were
% presented with the same set of stimuli. Hence, we use a single stimulus
% file (dataStim.mat) that applies to all subjects. This is compatible
% with scenarios with randomise presentation orders. In that case, the
% EEG/MEG trials should be sorted to match the single stimulus file. 
% The original order is preserved in a specific CND variable. If distinct
% subjects were presented with different stimuli, it is necessary to
% include a stimulus file per participant.
%
% CNSP-Workshop 2022
% https://cnspworkshop.net
% Author: Giovanni M. Di Liberto
% Copyright 2021 - Giovanni Di Liberto
%                  Nathaniel Zuk
%                  Michael Crosse
%                  Aaron Nidiffer
%                  Giorgia Cantisani
%                  (see license file for details)
% Last update: 24 June 2022

clear all
close all

addpath ../libs/cnsp_utils
addpath ../libs/cnsp_utils/cnd
addpath ../libs/mTRF-Toolbox_v2/mtrf
addpath ../libs/NoiseTools
addpath ../libs/eeglab
eeglab

datasets = {'../datasets/LalorNatSpeech/','../datasets/diliBach/','../datasets/cocktailAttSwitch/','../datasets/musicImagery/','../datasets/emilyDialogue/'};

for iiDataset = 1:length(datasets)
    iiDataset 

    % Parameters - Natural speech listening experiment
    dataMainFolder = datasets{iiDataset}; %'../datasets/LalorNatSpeech/';
    dataCNDSubfolder = 'dataCND/';
    
    eegFilenames = dir([dataMainFolder,dataCNDSubfolder,'dataSub*.mat']);
    nSubs = length(eegFilenames);
    
    % Average EEG data across participants
    for sub = 1:nSubs
        % Loading preprocessed EEG
        eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_',eegFilenames(sub).name];
        disp(['Loading preprocessed EEG data: ',eegFilenames(sub).name])
        load(eegPreFilename,'eeg')
        eeg = cndNormalise(eeg);
    
        if sub == 1
            eegAvg = eeg;
        else
            for iiTr = 1:length(eegAvg)
                % Sum EEG data
                minLen = min(size(eegAvg.data{iiTr},1),size(eeg.data{iiTr},1));
                eegAvg.data{iiTr} = eegAvg.data{iiTr}(1:minLen,:) + eeg.data{iiTr}(1:minLen,:);
                % Sum EEG external channels
                for iiExt = 1:length(eegAvg.extChan)
                    eegAvg.extChan{iiExt}.data{iiTr} = eegAvg.extChan{iiExt}.data{iiTr}(1:minLen,:) + eeg.extChan{iiExt}.data{iiTr}(1:minLen,:);
                end
            end
        end
    end
    
    eeg = eegAvg; clear eegAvg
    eeg = cndNormalise(eeg);
    
    eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_dataSubAvg.mat'];
    disp('Saving avg EEG file')
    save(eegPreFilename,'eeg')
    
    disp('Done')
end
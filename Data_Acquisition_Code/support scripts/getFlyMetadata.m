% getFlyMetadata.m
%
% Function that prompts user to enter metadata about the experimental fly
%  and returns it in a single struct. Also, saves identifying info about
%  fly (date and fly number).
% Metadata:
%   genotype - as string, fly genotype
%   manipulation - as string, anything done to fly (e.g. starvation,
%       pharmacology, etc.)
%   prepType - as string, how was fly mounted (which mount, postion, etc.)
%   dissectionNotes - as string, any comments about the dissection
%   ageUnits - whether age expressed in hours or days
%   age - in units specified by ageUnits, as single number or range 
%
% INPUTS:
%   dateDir - name of date directory, only YYMMDD portion, not full path
%   flyDir - name of fly directory, only fly## portion
%
% OUTPUTS:
%   flyData - struct of all of the above info (metadata and inputs)
%
% Created:  07/27/2018 - HY
%           02/02/2021 - MC
%           05/12/2021 - MC changed variables
%           08/30/2021 - MC added preselected flies
%           11/01/2021 - MC added sex field

function flyData = getFlyMetadata(dateDir, flyDir)
    flyData.dateDir = dateDir;
    flyData.flyDir = flyDir;
    listDir = 'C:\Users\wilson\Dropbox (HMS)\MATLAB\Ephys_Code\support scripts\';
    
    % load exhisting fly genotypes
    load([listDir 'myFlyList.mat'])
    % prompt user to select fly genotype
    [indx,~] = listdlg('ListString',myFlyList,'ListSize',[300, 500]);
    
    % pull genotype name
    genotype = myFlyList{indx};
    
    % add new entry for new genotypes
    if contains(genotype,'new','IgnoreCase',true)
        genotype = input('\n[INPUT] Enter full genotype: ', 's');
        myFlyList{length(myFlyList)+1} = genotype;
        
        %update list
        save([listDir 'myFlyList.mat'],'myFlyList')
    end
    
    
    % prompt user for manually entered metadata, in dialog box
    prompt = {'Sex','Light Cycle','Age (days since eclosed)','Starvation (hours)','Manipulations'};
    title = 'Fly Metadata';
    % dimensions; allows multiple lines and arbitrary length for genotype, 
    %  manipulation, prep type, dissection notes fields
    dims = [1 60; 1 60; 1 60; 1 60; 1 60];
    
    % do it this way so cancel just brings up dialog box again
    metadata = {};
    while isempty(metadata)
        metadata = inputdlg(prompt, title, dims);
    end
    
    % save metadata into flyData struct
    flyData.genotype = genotype;
    flyData.sex = metadata{1};
    flyData.lightcycle = metadata{2};
    flyData.age = str2double(metadata{3}); % convert age to numerical value
    flyData.starvation = metadata{4};
    flyData.notes = metadata{5};
    
end

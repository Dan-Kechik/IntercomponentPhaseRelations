source = fullfile(pwd, 'In'); %Datasets folder
%Get files from all datasets.
files = dir(fullfile(source, 'MPZ', '*', '*', '*.mat')); %'MPZ dataset'
files = dir(fullfile(source, '*', '*.mat'));
addPathes()

for ai = 1:numel(files)
    fullName = fullfile(files(ai).folder, files(ai).name);
    load(fullName, 'myFile');
    %Get phase data for the current file and update structure.
    myFile = processFileStruct(myFile);
    %Get conditions if they are - check the current datatset folder.
    dtstRoot = strsplit(strrep(files(ai).folder, source, ''), filesep);
    dtstRoot = fullfile(source, dtstRoot{2});
    if exist(fullfile(dtstRoot, 'conditions.m'), 'file')
        currentFolder = pwd;
        cd(dtstRoot)
        myFile.conditions = conditions(fullName);
        cd(currentFolder);
    else
        myFile.conditions = [];
    end
    %Update structure on disk.
    save(fullName, 'myFile', '-v7.3');
end
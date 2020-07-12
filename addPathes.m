function addPathes()
    if ~exist('ICPRcombinations', 'file') || ~exist('polyHarmonic', 'file')
        Root = fileparts(mfilename('fullpath'));
        gp = strsplit(genpath(Root), ';')';
        idxsGit = cellfun(@(x) fildInPath(x, '.git'), gp);
        idxsIn = cellfun(@(x) fildInPath(x, 'In'), gp);
        idxsOut = cellfun(@(x) fildInPath(x, 'Out'), gp);
        idxs = idxsGit | idxsIn | idxsOut;
        addpath(strjoin(gp(~idxs), ';')); 
    end
end

function pos = fildInPath(path, folder)
    dirs = strsplit(path, filesep);
    pos = nnz(strcmp(dirs, folder));
    pos = logical(pos);
end
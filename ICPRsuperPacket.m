%ICPRsuperPacket

%Put here dataset names for processing (from In folder).
dataset = {'Misalignment_vert_parall', 'Misalignment_norm', 'Misalignment_horz_parall'};
baseFreqz = {'shaftFreq'};
%Sign ICPRs for pictures saving.
icprLabelz = {' - PQI12', ' - PQI13', ' - PhI123'};
%Put desirable domains here.
domainz = {'acceleration', 'velocity', 'displacement'};
%See myFile.(domainz).(baseFreq).harmonicsIndexes to estimate index of desired ICPRs.
ICPRindexez = [1, 2, 7];

for cci = 1:numel(domainz)
    domain = domainz{cci};
    for aai = 1:numel(baseFreqz)
        baseFreq = baseFreqz{aai};
        for bbi = 1:numel(ICPRindexez)
            icprLabel = icprLabelz{bbi};
            ICPRindex = ICPRindexez(bbi);
            plotICPRstatistics
        end
        for ddi = 1:numel(dataset)
            currentFiles = strcmp(dataset{ddi}, {files.dataset});
            plotICPRsAndSpec(files(currentFiles), domain, baseFreq, icprLabelz, ICPRindexez);
        end
    end
end

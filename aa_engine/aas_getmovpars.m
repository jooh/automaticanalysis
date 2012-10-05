% return the combined movement parameters in mm (columns 1:3) and degrees
% (4:6) for a set of sessions (indices). Useful for re-plotting.
% movpar = aas_getmovpars(aap,sub,sessions)
function movpar = aas_getmovpars(aap,sub,sessions)

for s = 1:length(sessions)
    if s == 1
        P = spm_vol(aas_getimages(aap,sub,sessions(s),''));
    else
        P = [P; spm_vol(aas_getimages(aap,sub,sessions(s),''))];
    end
end
Params = zeros(numel(P),12);
for x = 1:numel(P)
    Params(x,:) = spm_imatrix(P(x).mat/P(1).mat);
end
movpar = [Params(:,1:3) Params(:,4:6)*180/pi];

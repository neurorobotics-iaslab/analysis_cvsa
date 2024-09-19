function [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, fix_event, cue_events, cf_event)
disp('   [INFO] extracting data infoo');
fixPOS = header.EVENT.POS(header.EVENT.TYP == fix_event);
fixDUR = header.EVENT.DUR(header.EVENT.TYP == fix_event);

cfPOS = header.EVENT.POS(header.EVENT.TYP == cf_event);
cfDUR = header.EVENT.DUR(header.EVENT.TYP == cf_event);

cuePOS = header.EVENT.POS(ismember(header.EVENT.TYP, cue_events));
cueDUR = header.EVENT.DUR(ismember(header.EVENT.TYP, cue_events));
cueTYP = header.EVENT.TYP(ismember(header.EVENT.TYP, cue_events));

cuePOS = cuePOS(length(cuePOS)-length(cfPOS)+1:end); % do this beacuse some value are in the eye calibration (we know they are at the beginning)
cueDUR = cueDUR(length(cueDUR)-length(cfDUR)+1:end);
cueTYP = cueTYP(length(cueTYP)-length(cfDUR)+1:end);

n_trial = sum(header.EVENT.TYP == 781);
end
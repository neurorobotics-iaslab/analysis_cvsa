function log_step(fmt, varargin)
% LOG_STEP  Prefixed printf used by every stage of the simulation so the
%           console transcript reads as a pipeline log.
    fprintf('[sim] ');
    fprintf(fmt, varargin{:});
    if isempty(fmt) || fmt(end) ~= newline
        fprintf('\n');
    end
end

function [proc_aligned, art_aligned, n_valid] = align_streams(proc, art, bufsize_proc, bufsize_art, chunk_size)
% ALIGN_STREAMS  The artifact ringbuffer (bufsize_art samples) fills before
%   the processing ringbuffer (bufsize_proc samples). Both streams are
%   indexed by chunk, but their "first valid chunk" differs by
%   (bufsize_proc - bufsize_art) / chunk_size. This trims both streams to
%   start at the chunk where BOTH are valid (= processing's first valid
%   chunk) so they can be consumed in lockstep.
%
%   Inputs
%       proc           [n_chunks x n_features]  (NaN before processing valid)
%       art            [n_chunks x 1] logical    (false before artifact valid)
%       bufsize_proc   samples
%       bufsize_art    samples
%       chunk_size     samples per chunk
%
%   Outputs
%       proc_aligned, art_aligned : trimmed streams starting from first
%                                   chunk where both are valid
%       n_valid : number of chunks in the aligned streams
    first_valid_proc = ceil(bufsize_proc / chunk_size);   % 1-indexed
    proc_aligned     = proc(first_valid_proc:end, :);
    art_aligned      = art (first_valid_proc:end);
    n_valid          = size(proc_aligned, 1);
end

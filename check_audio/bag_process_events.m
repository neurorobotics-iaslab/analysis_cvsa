function [events_out] = bag_process_events(bag)
    events = bag.events;
    events.TYP = events.TYP;
    events.time = events.time;
    event_types = unique(events.TYP(events.TYP<0x8000))';
    events_out.TYP = [];
    events_out.time = [];
    events_out.dur = [];
    for ev = event_types
        start_idx = (events.TYP == ev);
        end_idx = (events.TYP==ev+0x8000);
        events_out.TYP = [events_out.TYP; events.TYP(start_idx)];        
        events_out.dur = [events_out.dur; (events.time(end_idx) - events.time(start_idx))];
        events_out.time = [events_out.time; events.time(start_idx)];
    end
    [events_out.time,idx] = sort(events_out.time);
    events_out.TYP = events_out.TYP(idx);
    events_out.dur = events_out.dur(idx);
    events_out.POS = arrayfun(@(x) find(abs(bag.audio.time - x) == min(abs(bag.audio.time  - x)), 1), events_out.time);
    events_out.DUR = floor(events_out.dur*bag.audio.fs);
end
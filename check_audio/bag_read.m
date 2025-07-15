function [out] = bag_read(filename)
    bag = rosbag(filename);
    audio_bag = bag.select('Topic','/audio');
    events_bag = bag.select('MessageType','rosneuro_msgs/NeuroEvent','Topic','/events/bus');
 
    audio_msgs = audio_bag.readMessages('DataFormat','struct');
    events_msgs = events_bag.readMessages('DataFormat','struct');

    msg_size = size(audio_msgs{1}.Audio.Data,1);

    data = zeros(size(audio_msgs,1)*msg_size,1);
    audio_stamps = zeros(size(audio_msgs,1),1);
    
    for i = 1:size(audio_msgs,1)
        data(1+(i-1)*msg_size:(i)*msg_size,:) = audio_msgs{i}.Audio.Data;
        audio_stamps(i) = audio_bag.MessageList.Time(i);
    end
    
    events_codes = zeros(size(events_msgs,1),1);
    for i = 1:size(events_msgs,1)
        events_codes(i) = events_msgs{i}.Event;
    end
    
    audio.data = data(msg_size:end);
    audio.data = resample(audio.data,44112,cast(audio_msgs{1}.SampleRate,'double')); % resampling is necessary for the very first sessions since they were originally recorded at a different sampling rate has no effect of the others
    audio.data = int16(audio.data);
    audio.stamps = audio_stamps;
    audio.time = linspace(0,audio_bag.MessageList.Time(end)-audio_bag.MessageList.Time(1),size(audio.data,1));
    audio.fs = 44112;
    audio.bit_depth = cast(16,'double');
    audio.chunk_size = cast(msg_size,'double');

    events.TYP = events_codes;
    events.time = events_bag.MessageList.Time -audio_bag.MessageList.Time(1);
    
    out.audio = audio;
    out.events = events;
end
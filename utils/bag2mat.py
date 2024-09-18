import rosbag
import rospy
from scipy.io import savemat
import os
import glob

class stamp:
    def __init__(self, secs, nsecs):
        self.secs = secs
        self.nsecs = nsecs

class header:
    def __init__(self, seq, stamp, frame_id):
        self.seq = seq
        self.nsecs = stamp
        self.frame_id = frame_id

class neuroheader:
    def __init__(self, seq):
        self.seq = seq

class left_pupil:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
class right_pupil:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class message:
    def __init__(self, header, left_pupil, right_pupil):
        self.header = header
        self.left_pupil = left_pupil
        self.right_pupil = right_pupil

class neuroEvent:
    def __init__(self, header, neuroheader, event, family, description, duration, version):
        self.header = header
        self.neuroheader = neuroheader
        self.event = event
        self.family = family
        self.description = description
        self.duration = duration
        self.version = version


def convert_bag_to_mat(bag_file, output_mat_file):
    bag = rosbag.Bag(bag_file)
    messages = {'cvsa_pupils': [], 'events_bus' : []}

    for topic, msg, t in bag.read_messages():
        # Convert ROS messages to dictionary
        if topic == '/cvsa/pupils':
            #cstamp = stamp(msg.header.stamp.secs, msg.header.stamp.nsecs)
            cstamp = msg.header.stamp.secs*1000000000 + msg.header.stamp.nsecs
            cheader = header(msg.header.seq, cstamp, msg.header.frame_id)
            cleft_pupil = left_pupil(msg.left_pupil.x, msg.left_pupil.y, msg.left_pupil.z)
            cright_pupil = right_pupil(msg.right_pupil.x, msg.right_pupil.y, msg.right_pupil.z)
            messages['cvsa_pupils'].append(message(cheader, cleft_pupil, cright_pupil))
        if topic == '/events/bus':
            #cstamp = stamp(msg.header.stamp.secs, msg.header.stamp.nsecs)
            cstamp = msg.header.stamp.secs*1000000000 + msg.header.stamp.nsecs
            cheader = header(msg.header.seq, cstamp, msg.header.frame_id)
            cneuroheader = neuroheader(msg.neuroheader.seq)
            event = msg.event
            family = msg.family
            description = msg.description
            duration = msg.duration
            version = msg.version
            messages['events_bus'].append(neuroEvent(cheader, cneuroheader, event, family, description, duration, version))



    bag.close()
    print('Saving file: '+ output_mat_file)
    savemat(output_mat_file, messages)

    

if __name__ == "__main__":
    directory = '/home/paolo/cvsa_ws/record/bag'
    bag_files = glob.glob(os.path.join(directory, '*.bag'))

    for bag_file in bag_files:
        print('Loading file: ' + bag_file)
        output_mat_file = bag_file[0:len(bag_file)-3] + 'mat'  
        convert_bag_to_mat(bag_file, output_mat_file)

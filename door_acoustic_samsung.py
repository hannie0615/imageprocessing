import numpy as np

class door_acoustic_samsung:
    list_arg_max_freq = np.zeros(40)

    door_open_ground_truth = np.array(
        [65,65,65,65,1000,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,1000,132,132,132])
    door_close_ground_truth = np.array(
        [132,132,132,132,132,132,132,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,97,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,65,65])
    door_error_ground_truth = np.array(
        [97,97,97,97,97,97,97,82,82,82,82,82,82,82,82,82,82,97,97,97,97,97,97,97,97,97,82,82,82,82,82,82,82,82,82,82])
    door_error_big_ground_truth = np.array(
        [65,65,65,65,65,65,65,65,65,65,65,65,65,65,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,65,65,65,65,65,65,65,65])
    attack_ground_truth = np.array([86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,99,99])
    def __init__(self):
        self.list_arg_max_freq=np.zeros(56)
    def find_max_freq(self,input_data):
        self.list_arg_max_freq[0:40] = self.list_arg_max_freq[16:56]
        for ii in range(16):
            temp_data=input_data[ii*256:ii*256+512]
            fft_out = np.fft.fft(temp_data,512)
            abs_out = abs(fft_out)
            self.list_arg_max_freq[40+ii] = np.argmax(abs_out[35:250]) + 35



    def check_door_open(self,index):
        sum_of_all_difference=0
        #for ii in range(len(self.door_open_ground_truth)):
        for ii in range(40):
            if self.door_open_ground_truth[ii]==1000:
                temp_diff=0
            elif self.door_open_ground_truth[ii] == 65 and self.list_arg_max_freq[ii+index] == 129:
                temp_diff = 0
            else:
                temp_diff=abs(self.door_open_ground_truth[ii]-self.list_arg_max_freq[ii+index])
            sum_of_all_difference+=temp_diff

        #print(sum_of_all_difference)
        if sum_of_all_difference<200:
            return 60

    def check_door_error(self,index):
        sum_of_all_difference = 0

        sum_of_all_difference = sum(abs(self.door_error_ground_truth - self.list_arg_max_freq[4+index:40+index]))
        #print(sum_of_all_difference)
        if sum_of_all_difference < 50:
            return 70

    def check_door_error_big(self,index):
        sum_of_all_difference = 0
        #for ii in range(len(self.door_error_big_ground_truth)):
        for ii in range(36):
            if self.door_error_big_ground_truth[ii] == 65 and self.list_arg_max_freq[ii+index] == 129:
                temp_diff = 0
            else:
                temp_diff = abs(self.door_error_big_ground_truth[ii] - self.list_arg_max_freq[ii+index])
            sum_of_all_difference += temp_diff
        #print(sum_of_all_difference)
        if sum_of_all_difference < 150:
            return 80



    def check_door_close(self,index):
        sum_of_all_difference = 0
        #for ii in range(len(self.door_close_ground_truth)):
        """for ii in range(40):
            temp_diff = abs(self.door_close_ground_truth[ii] - self.list_arg_max_freq[ii])
            sum_of_all_difference += temp_diff
        """
        sum_of_all_difference = sum(abs(self.door_close_ground_truth - self.list_arg_max_freq[index:index+40]))
        #print(sum_of_all_difference)
        if sum_of_all_difference < 150:
            return 90


    def check_attack(self,index):
        sum_of_all_difference = 0
        #for ii in range(len(self.attack_ground_truth)):
        """for ii in range(29):
            temp_diff = abs(self.attack_ground_truth[ii] - self.list_arg_max_freq[ii + 11])
            sum_of_all_difference += temp_diff
        """
        sum_of_all_difference = sum(abs(self.attack_ground_truth - self.list_arg_max_freq[11+index:40+index]))
        #print(sum_of_all_difference)
        if sum_of_all_difference < 50:
            return 100


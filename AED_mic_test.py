
import pyaudio
from six.moves import queue
import time
import numpy as np
import pickle
from scipy import signal
import python_speech_features
import door_acoustic_samsung
import AED_CNN_6

# 녹음용 값
# 16khz
RATE = 16000
# 버퍼
CHUNK = 4096



class AED_CNN_7:

    def __init__(self):
        #f=open('weight_data_1029_1.bin', 'rb')
        #f=open('weight_data_1030_1.bin','rb')
        f = open('weight_data_1109_2.bin', 'rb')
        self.w1 = pickle.load(f)
        self.w2 = pickle.load(f)
        self.w3 = pickle.load(f)
        self.w4 = pickle.load(f)
        self.w5 = pickle.load(f)
        self.b1 = pickle.load(f)
        self.b2 = pickle.load(f)
        self.b3 = pickle.load(f)
        self.b4 = pickle.load(f)
        self.b5 = pickle.load(f)
        f.close()
        
        #f = open('weight_01_1104_1.bin', 'rb')
        f = open('weight_01_1109_2.bin', 'rb')
        #f = open('weight_01_1029_2_1.bin','rb')
        self.w1_01 = pickle.load(f)
        self.w2_01 = pickle.load(f)
        self.w3_01 = pickle.load(f)
        self.w4_01 = pickle.load(f)
        self.b1_01 = pickle.load(f)
        self.b2_01 = pickle.load(f)
        self.b3_01 = pickle.load(f)
        self.b4_01 = pickle.load(f)
        f.close()


    # def ToMono(self, y_stereo):
    #     return (y_stereo[0] + y_stereo[1]) / 2


    def find_01(self, fea_output):

        L1_01 = np.zeros((8, 8, 4))
        for ii in range(4):
            tmp_a1 = self.w1_01[:, :, 0, ii]
            tmp_a1.reshape(3, 3)
            tmp_a1 = np.flip(tmp_a1)
            tmp_filtered_output = signal.convolve2d(fea_output, tmp_a1, mode='same')
            tmp_filtered_output = tmp_filtered_output + self.b1_01[ii]
            tmp_filtered_output = np.maximum(tmp_filtered_output, 0)

            tmp_max_pooling = np.zeros((8, 8))
            
            for kk in range(8):
                for jj in range(8):
                    tmp_max_pooling[kk, jj] = np.max(tmp_filtered_output[kk * 4:kk * 4 + 4, jj * 4:jj * 4 + 4])
            L1_01[:, :, ii] = tmp_max_pooling



        L2_01 = np.zeros((2, 2, 8))
        for ii in range(8):
            tmp_sum_filtered_output = np.zeros((8, 8))
            for jj in range(4):
                tmp_L1 = L1_01[:, :, jj]
                tmp_L1 = tmp_L1.reshape(8, 8)
                tmp_a2 = self.w2_01[:, :, jj, ii]
                tmp_a2.reshape(3, 3)
                tmp_a2 = np.flip(tmp_a2)
                tmp_filtered_output = signal.convolve2d(tmp_L1, tmp_a2, mode='same')
                tmp_sum_filtered_output = tmp_sum_filtered_output + tmp_filtered_output
                
            tmp_sum_filtered_output = tmp_sum_filtered_output + self.b2_01[ii]
            filtered_output = np.maximum(tmp_sum_filtered_output, 0)
            
            tmp_max_pooling = np.zeros((2, 2))

            for kk in range(2):
                for oo in range(2):
                    tmp_max_pooling[kk, oo] = np.max(filtered_output[kk * 4:kk * 4 + 4, oo * 4:oo * 4 + 4])
            L2_01[:, :, ii] = tmp_max_pooling
        L2_reshape = L2_01.reshape(1, 32)
        L3_01 = np.matmul(L2_reshape, self.w3_01)
        L3_01 = L3_01 + self.b3_01
        L3_01 = np.maximum(L3_01, 0)
        L4_01 = np.matmul(L3_01, self.w4_01)
        L4_01 = L4_01 + self.b4_01

        if L4_01<0:
            out=0
        else:
            out=1
        #print(L4_01)

        return out


    def find_event(self, y_mono, Door_test, rate=16000):
        # y_mono = self.ToMono(y_stereo)
        Door_test.find_max_freq(y_mono)
        for ii in range(16):
            if Door_test.check_door_open(ii) == 60:
                return 60,1
            if Door_test.check_door_error(ii) == 70:
                return 70,1
            if Door_test.check_door_error_big(ii) == 80:
                return 80,1
            if Door_test.check_door_close(ii) == 90:
                return 90,1
            if Door_test.check_attack(ii) == 100:
                return 100,1

        squred_sig = y_mono ** 2
        mean_squared = (np.mean(squred_sig))
        #print(mean_squared)
        
        if mean_squared < 0.000001:
            return 0, 3
        if max(abs(y_mono))!= 0:
            y_mono = y_mono / (max(abs(y_mono)) * 1.2)

        fea_output = python_speech_features.base.logfbank(y_mono, samplerate=16000, winlen=0.032, winstep=0.016,
                                                          nfilt=32, lowfreq =30, highfreq=4000)

        ooo=self.find_01(fea_output)

        # print(ooo)
        if ooo==0:
            return 0,6




        L1_out = np.zeros((16, 16, 4))
        for ii in range(4):
            tmp_a1 = self.w1[:, :, 0, ii]
            tmp_a1.reshape(3, 3)
            tmp_a1 = np.flip(tmp_a1)
            tmp_filtered_output = signal.convolve2d(fea_output, tmp_a1, mode='same')
            tmp_filtered_output = tmp_filtered_output + self.b1[ii]
            tmp_filtered_output = np.maximum(tmp_filtered_output, 0)

            tmp_max_pooling = np.zeros((16, 16))
            for kk in range(16):
                for jj in range(16):
                    tmp_max_pooling[kk, jj] = np.max(tmp_filtered_output[kk * 2:kk * 2 + 2, jj * 2:jj * 2 + 2])
            L1_out[:, :, ii] = tmp_max_pooling


        L2_out = np.zeros((8, 8, 8))
        for ii in range(8):
            tmp_sum_filtered_output = np.zeros((16, 16))
            for jj in range(4):
                tmp_L1 = L1_out[:, :, jj]
                tmp_L1 = tmp_L1.reshape(16, 16)
                tmp_a2 = self.w2[:, :, jj, ii]
                tmp_a2.reshape(3, 3)
                tmp_a2 = np.flip(tmp_a2)
                tmp_filtered_output = signal.convolve2d(tmp_L1, tmp_a2, mode='same')
                tmp_sum_filtered_output = tmp_sum_filtered_output + tmp_filtered_output
            tmp_sum_filtered_output = tmp_sum_filtered_output + self.b2[ii]
            filtered_output = np.maximum(tmp_sum_filtered_output, 0)

            tmp_max_pooling = np.zeros((8, 8))

            for kk in range(8):
                for oo in range(8):
                    tmp_max_pooling[kk, oo] = np.max(filtered_output[kk * 2:kk * 2 + 2, oo * 2:oo * 2 + 2])
            L2_out[:, :, ii] = tmp_max_pooling


        L3_out = np.zeros((2, 2, 16))
        for ii in range(16):
            tmp_sum_filtered_output = np.zeros((8, 8))
            for jj in range(8):
                tmp_L2 = L2_out[:, :, jj]
                tmp_L2 = tmp_L2.reshape(8, 8)
                tmp_a3 = self.w3[:, :, jj, ii]
                tmp_a3.reshape(3, 3)
                tmp_a3 = np.flip(tmp_a3)
                tmp_filtered_output = signal.convolve2d(tmp_L2, tmp_a3, mode='same')
                tmp_sum_filtered_output = tmp_sum_filtered_output + tmp_filtered_output
            tmp_sum_filtered_output = tmp_sum_filtered_output + self.b3[ii]
            filtered_output = np.maximum(tmp_sum_filtered_output, 0)
            tmp_max_pooling = np.zeros((2, 2))
            for kk in range(2):
                for oo in range(2):
                    tmp_max_pooling[kk, oo] = np.max(filtered_output[kk * 4:kk * 4 + 4, oo * 4:oo * 4 + 4])
            L3_out[:, :, ii] = tmp_max_pooling



        L3_reshape = L3_out.reshape(1, 64)
        L4_out = np.matmul(L3_reshape, self.w4)
        L4_out = L4_out + self.b4
        L4_out = np.maximum(L4_out, 0)
        L5_out = np.matmul(L4_out, self.w5)
        L5_out = L5_out + self.b5
        out=np.argmax(L5_out)

        tmp_softmax=np.exp(L5_out)

        softmax_out=tmp_softmax[0,np.argmax(L5_out)]/np.sum(tmp_softmax)


        return out, softmax_out



class MicrophoneStream(object):
    """마이크 입력 클래스"""

    def __init__(self, rate, chunk):
        self._rate = rate
        self._chunk = chunk

        # 마이크 입력 버퍼 생성
        self._buff = queue.Queue()
        self.closed = True

    # 클래스 열면 발생함.
    def __enter__(self):
        # pyaudio 인터페이스 생성
        self._audio_interface = pyaudio.PyAudio()
        # 16비트, 모노로 마이크 열기
        # 여기서 _fill_buffer 함수가 바로 callback함수 인데
        # 실제 버퍼가 쌓이면 이곳이 호출된다.
        # 즉, _fill_buffer 마이크 입력을 _fill_buffer 콜백함수로 전달 받음
        self._audio_stream = self._audio_interface.open(
            format=pyaudio.paInt16,
            channels=1, rate=self._rate,
            input=True, frames_per_buffer=self._chunk,
            stream_callback=self._fill_buffer,
        )
        self.closed = False
        return self

    def __exit__(self, type, value, traceback):
        # 클래스 종료시 발생
        # pyaudio 종료
        self._audio_stream.stop_stream()
        self._audio_stream.close()

        self.closed = True
        # Signal the generator to terminate so that the client's
        # streaming_recognize method will not block the process termination.
        self._buff.put(None)
        self._audio_interface.terminate()

    # 마이크 버퍼가 쌓이면(CHUNK = 1600) 이 함수 호출 됨.
    def _fill_buffer(self, in_data, frame_count, time_info, status_flags):
        # 마이크 입력 받으면 큐에 넣고 리턴
        self._buff.put(in_data)
        return None, pyaudio.paContinue

    # 제너레이터 함수
    def generator(self):
        # 클래스 종료될 떄까지 무한 루프 돌림
        while not self.closed:

            # 큐에 데이터를 기다림.
            # block 상태임.
            chunk = self._buff.get()

            # 데이터가 없다면 문제 있음
            if chunk is None:
                return

            # data에 마이크 입력 받기
            data = [chunk]

            # 추가로 받을 마이크 데이터가 있는지 체크
            while True:
                try:
                    # 데이터가 더 있는지 체크
                    chunk = self._buff.get(block=False)
                    if chunk is None:
                        return
                    # 데이터 추가
                    data.append(chunk)
                except queue.Empty:
                    # 큐에 데이터가 더이상 없다면 break
                    break

            # 마이크 데이터를 리턴해줌
            yield b''.join(data)
# [END audio_stream]


def main():
    Door = door_acoustic_samsung.door_acoustic_samsung()
    AED_CNN = AED_CNN_7()
    outlist = [0, 0, 0, 0]
    temp_data = np.zeros(8192+256)
    
    # 마이크 열기
    with MicrophoneStream(RATE, CHUNK) as stream:
        # 마이크 데이터 핸들을 가져옴
        audio_generator = stream.generator()
        
        for i in range(100): # 1000번

            print(audio_generator)

            for x in audio_generator: # x = 8192
                
                # 마이크 음성 데이터
                temp_data[0:4096+256] = temp_data[4096:4096+4096+256]
                
                for jj in range(4096):
                    temp_data[4096+256+jj]=int.from_bytes(x[jj*2:jj*2+2], byteorder='big', signed=True)/32768
                
                ## AED_CNN 적용 
                out, out2 = AED_CNN.find_event(temp_data, Door, 16000)
                print(out, out2)
                
                if out > 50:
                    print(out)
                else:
                    outlist.append(out)
                    del outlist[0]

                final_out = 0
                temp_out = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                
                for ii in range(3):
                    temp_out[outlist[ii]] += 1

                for ii in range(10):
                    if temp_out[ii] > 1.9:
                        final_out = ii
                        break

                print(final_out)
                break



if __name__ == '__main__':
    main()
    
    
    
    
    
    
    
    
    
    
    
    
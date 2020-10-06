import numpy as np
import matplotlib.pyplot as plt
import sys
from read_geometries import read_geoms
from tqdm import tqdm
from tidynamics import acf
from tidynamics import correlation

class Correlation_Functions:
    """
    A class for calculating correlation functions from MD simulation data. Convenience
    functions for reading in xyz files are provided, but all methods can be called as long as
    sequences of arrays are provided. Since virtually all correlation functions from MD will
    be between arrays of 3D vectors, I provide this as the default.
    """
    def __init__(self):
        self.signal = None
        self.frequencies = None
        self.frame_spacing = None

    def acf(self, A: np.ndarray, frame_spacing: float):
        """
        Computes the autocorrelation of data which is provided as a 2D numpy array. Since the 
        components are all independent, this is flatted into a sequence of 1D arrays, and the 
        correlation function is computed for each component of each atom and averaged.
        The correlation function is computed by taking the FFT and inverse FFT fo the data.

        Takes the spacing in time between frames in femtoseconds.
        """
        self.frame_spacing = frame_spacing
        data = np.asarray([a.ravel() for a in A])
        data = np.vstack(data.T)
        self.signal = np.zeros(data.shape[1] // 2)
        N = data.shape[1]
        M = len(self.signal)
        for i in tqdm(range(N - M)):
            for iAtom in range(data.shape[0]):
                self.signal += acf(data[iAtom, i:M+i].T)
        self.signal /= (data.shape[0])
        self.signal /= (self.signal[0])
        self.frequencies = np.fft.fftfreq(self.signal.size, self.frame_spacing)[0:int(self.signal.size / 2)]
        self.frequencies *= 10**(15) / (2.9979*(10**10))
        return self.frequencies, self.signal
    
    def block_average_acf(self, A: np.ndarray, frame_spacing: float, block_width: int):
        """
        Computes autocorrelation over as many blocks of size width as possible.
        Simply calls acf and averages the results over all blocks.
        """
        num_blocks = int(len(A) / block_width) - 1
        accumulated_signal = np.zeros(int(block_width / 2))
        for i in range(num_blocks):
            _, signal = self.acf(A[i*block_width:(i+1)*block_width], frame_spacing)
            accumulated_signal += signal
        accumulated_signal /= num_blocks
        self.signal = accumulated_signal

    def plot_signal(self):
        plt.plot(np.asarray(range(len(self.signal))) * self.frame_spacing, self.signal)
        plt.xlabel("Time (fs)")
        plt.ylabel("Signal Amplitude (a.u.)")
        plt.show()
    
    def plot_fourier_transformed_signal(self):
        vdos = np.fft.fft(self.signal).real
        plt.plot(self.frequencies, vdos[0:int(vdos.size / 2)])
        plt.xlabel("Energy ($cm^{-1}$)")
        plt.ylabel("Intensity (a.u.)")
        plt.ylim(0,15)
        plt.show()

if __name__ == '__main__':
    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get input file.")
        sys.exit(1)
    corr = Correlation_Functions()
    _, _, velocities = read_geoms(ifile)
    corr.block_average_acf(velocities, 2.5, 5000)
    corr.plot_signal()
    corr.plot_fourier_transformed_signal()


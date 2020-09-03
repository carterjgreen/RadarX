import numpy as np
from scipy.constants import c, pi

class Target:
    '''
    Base Class for a target that will be used with the ADC class
    
    Parameters
    -----------
    R0 : Initial range to target (m)
    vr0 : Target's initial closing radial velocity (m/s)
    vx0 : Target's cross range velocity (m/s)
    phi0 : Target's initial azimuth angle 
    EL : Target's elevation (1, 2, or 3)
    t_ap : Time at which the target appears (s)
    t_dis : Time at which the target disappears (s)
    CPI : Target's coherent processing interval (s)
    '''
    def __init__(self, R0, phi0, vr0, vx0=0, EL=1, t_ap=0, t_dis=30, CPI=5, SNR=10):
        self.R0 = R0
        self.phi0 = phi0
        self.R = R0
        self.phi = phi0
        self.vr = vr0
        self.vx = vx0
        self.EL = EL
        self.t_ap = t_ap
        self.t_dis = t_dis
        self.CPI = CPI
        self.SNR = SNR
        

class ADC:
    """
    ADC object
    
    Parameters
    ----------
    f0 : 
        RF center frequency (Hz)
        
    Rdmax : 
        Maximum detection range (m)
        
    noise_var : 
        variance of Gaussian I and Q voltage samples (V^2)
        
    tau : 
        Pulse width (s)
        
    PRF : 
        Pulse repetition frequency (Hz)
        
    PRI_counter : int
        number of pulses, so far, at current PRF
        
    theta0 : 
        target phase, needed for coherence 
    
    PSI : 
        Phase progression to apply in azimuth (degrees)
        
    EL : 
        Elevation beam (1, 2, or 3)
        
    sim_time : 
        Simulation time (s)
        
    """
    light = c
    
    def __init__(self, f0=10e9, Rmax=30e3, noise_var=0.05):
        self.f0 = f0
        self.lambd = ADC.light/self.f0
        self.k = 2*pi/self.lambd
        self.Rmax = Rmax
        self.noise_var = noise_var
        self.tau = 1e-6
        self.Rmin = ADC.light*self.tau/2
        self.PRF = 17e3
        self.PRI = 1/self.PRF
        self.PRI_counter = 0
        self.previous_PRF = 0
        self.theta0 = 0
        self.PSI = 0
        self.EL = 1
        self.time = 0
        self.lookAngle = 0
        self.L = 0
        
        self.targets = []
    
        
    def add_target(self, target):
        self.targets.append(target)
        
        
    def update_targets(self):
        for target in self.targets:
            if (target.t_ap <= self.time < target.t_dis):
                target.R = target.R0 - (self.time-target.t_ap)*target.vr
                target.phi += target.phi0 - np.arctan((self.time-target.t_ap)*np.divide(target.vx, target.R))
                
                if target.R < 0:
                    target.t_dis = self.time
                
                
    def visible_targets(self):
        # Return an array containing ranges of all visible targets
        return np.array([target.R for target in self.targets if (target.t_ap <= self.time < target.t_dis) and (target.EL == self.EL)])
    
    
    def visible_SNR(self):
        # Return an array containing the SNR of all visible targets
        return np.array([target.SNR for target in self.targets if target.t_ap <= self.time < target.t_dis and (target.EL == self.EL)])
    
    
    def visible_CPI(self):
        # Return an array containing the CPI of all visible targets
        return np.array([target.CPI for target in self.targets if target.t_ap <= self.time < target.t_dis and (target.EL == self.EL)])        
          
        
    def set_PRF(self, PRF_new):
        self.PRF = PRF_new
        Ru = ADC.light/(2*self.PRF)
        delta_Rmin = ADC.light*self.tau/2
        self.L = np.int(Ru/delta_Rmin)
        
        
    def range2bin(self):
        Ru = ADC.light/(2*self.PRF)
        delta_Rmin = ADC.light*self.tau/2
        visible = self.visible_targets()
        self.L = np.int(Ru/delta_Rmin)
        l_targets = np.floor(visible/ delta_Rmin)
        l_app = np.floor(np.remainder(l_targets, self.L)).astype('int')
        
        return l_app
        
        
    def iq_noise(self):
        iNoise = np.random.normal(0, np.sqrt(self.noise_var), (1, self.L))
        qNoise = np.random.normal(0, np.sqrt(self.noise_var), (1, self.L))
        
        return iNoise, qNoise
    
    
    def clutter(self, on=True, clutter_power=10):
        if on and self.EL == 1:
            clutter = np.ones((1, self.L)) * np.sqrt(np.power(clutter_power/10, 10)*self.noise_var*2)
        else:
            clutter = 0
        
        return clutter
            
    
    def check_coherence(self):
        CPI = self.visible_CPI()
        phase = self.theta0 * np.ones_like(CPI)
        ind = (self.PRI_counter-1)/self.PRF > CPI
        
        if ind.any():
            phase[ind] = np.random.uniform(-pi, pi, size=phase[ind].shape)
        
        return phase
        
        
    def create_return(self):
        target_bins = self.range2bin()
        bins = np.zeros((1, self.L), dtype='complex128')
        SNR = self.visible_SNR()
        d = self.visible_targets()
        phase = self.check_coherence()
        
        signal_mag = np.sqrt(np.power(10, SNR/10)*self.noise_var*2)
        signal_phase = np.remainder(-2*self.k*d + phase, 2*pi)
        signal = signal_mag * np.exp(1j*signal_phase)
        
        np.put(bins, target_bins, signal)
        (iNoise, qNoise) = self.iq_noise()
        clutter = self.clutter()
        
        I = np.real(bins) + iNoise + clutter
        Q = np.imag(bins) + qNoise

        I[0, 0], Q[0, 0] = 0, 0
        return (I, Q)
        
        
    def read(self):
        # Check PRI
        if self.previous_PRF != self.PRF:
            self.theta0 = float(np.random.uniform(-pi, pi, 1))
            self.PRI_counter = 0
            
        self.previous_PRF = self.PRF
        self.PRI_counter += 1
        self.time += self.PRI
        self.update_targets()
        return self.create_return()
        
        
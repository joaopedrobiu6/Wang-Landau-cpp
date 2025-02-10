# Metropolis algorithm simulation example for 2D Ising model.
#
# Copyright (c) 2022-2024 Jack Lidmar
# All rights reserved.
#

import numpy as np
from numpy.random import rand as rnd

from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Slider, TextBox

class Ising:

    Tc = 2.0 / np.log(1.0 + np.sqrt(2.0))

    def __init__(self,L,T):
        self.init(L,T)

    def init(self,L,T):
        self.L = L
        self.s = np.ones((L,L), dtype=int)   # All spin up
        self.M = L*L        # Magnetization
        self.E = - 2*L*L    # Energy
        self.draw = False
        self.reset(T)

    def reset(self,T):
        self.T = T
        self.Nsamp = 0
        self.t = 0
        self.accratio = 0
        self.m = 0
        self.absm = 0
        self.e = 0
        self.chi = 0
        self.cv = 0
        self.binder = 0
        self.dEmax = 12
        self.wtab = np.exp(-np.arange(-self.dEmax,+self.dEmax+1) / self.T)
    
    def mcmc(self,n):
        self.metropolis(n)

    def metropolis(self,n):
        "MCMC moves"
        acc = 0
        n *= self.L ** 2
        for ii in range(n):
            x = int(self.L * rnd())
            y = int(self.L * rnd())
            sxy = self.s[x,y]
            dE = 2*sxy*( 
                self.s[(x + 1) % self.L,y] +
                self.s[(x - 1 + self.L) % self.L,y] +
                self.s[x, (y + 1) % self.L] +
                self.s[x, (y - 1 + self.L) % self.L] )
            
            # Metropolis:
            # if dE <= 0 or np.exp(-dE/self.T) > rnd():
            if dE <= 0 or self.wtab[dE+self.dEmax] > rnd(): # Use a table
                self.s[x,y] = -sxy  # Accept
                self.E += dE
                self.M -= 2*sxy
                acc += 1
        self.accratio += acc/n

    def energy_and_magnetization(self):   # Calculate energy and magnetization.
        xm = self.L - 1
        E = 0
        M = 0
        for x in range(self.L):
            ym = self.L -1
            for y in range(self.L):
                E -= self.s[x,y]*(self.s[xm,y] + self.s[x,ym])
                M += self.s[x,y]
                ym = y
            xm = x
        return E, M


    "Wolff cluster update."
    def Wolff(self,n):
        self.numflips = 0
        for ii in range(n):
            x = int(self.L * rnd()) # Choose a random seed for
            y = int(self.L * rnd()) # the Wolff cluster
            self.Wolff_flip(x,y)

    def Wolff_flip(self,x,y):
        dM = 2*self.s[x,y]
        self.s[x,y] = -self.s[x,y]  # Unconditionaly flip the spin.
        self.numflips += 1  # Keep track of number of spins flipped.
        self.M -= dM        # and magnetization

        xp = (x + 1) % self.L           # right neigbour
        xm = (x - 1 + self.L) % self.L  # left neighbour
        yp = (y + 1) % self.L           # up
        ym = (y - 1 + self.L) % self.L  # down

        # Update the energy:
        self.E += dM*(self.s[xm, y] + self.s[xp, y] + self.s[x, yp] + self.s[x, ym])

        dE = dM*self.s[xp,y]
        if dE > 0 and rnd() < 1 - np.exp(-dE/self.T):
            self.Wolff_flip(xp,y)

        dE = dM*self.s[xm,y]
        if dE > 0 and rnd() < 1 - np.exp(-dE/self.T):
            self.Wolff_flip(xm,y)

        dE = dM*self.s[x,yp]
        if dE > 0 and rnd() < 1 - np.exp(-dE/self.T):
            self.Wolff_flip(x,yp)

        dE = dM*self.s[x,ym]
        if dE > 0 and rnd() < 1 - np.exp(-dE/self.T):
            self.Wolff_flip(x,ym)

    def sample(self):
        self.m += self.M
        self.absm += abs(self.M)
        self.e += self.E
        self.chi += self.M ** 2
        self.cv += self.E ** 2
        self.binder += self.M ** 4
        self.Nsamp += 1

    def run(self,nequil,nsweeps,sweepsPerSample):
        self.mcmc(nequil)
        for t in range(nsweeps):
            self.mcmc(sweepsPerSample)
            self.sample()
        
        # Check that there were no mistakes in energy calculation:
        E, M = self.energy_and_magnetization()
        if self.E != E:
            print(f"Energy missmatch: {self.E=}, {E=}")

        self.results()

    def results(self):
        self.m /= self.Nsamp
        self.absm /= self.Nsamp
        self.e /= self.Nsamp
        self.chi /= self.Nsamp
        self.cv /= self.Nsamp
        self.binder /= self.Nsamp
        self.binder = 1 - self.binder/(3 * self.chi**2)
        self.chi -= self.absm ** 2
        self.cv -= self.e ** 2

        self.m /= self.L ** 2
        self.absm /= self.L ** 2
        self.e /= self.L ** 2
        self.chi /= self.L ** 2 * self.T
        self.cv /= self.L ** 2 * self.T ** 2

    def write(self,file):
        file.write(f"{self.T} {self.m} {self.absm} {self.chi} {self.e} {self.cv} {self.binder} {self.L}\n")
    
    @staticmethod
    def animate(L : int, Ts : np.iterable, Nit = 1000000, Nskip = 10, Wolff = False):
        system = Ising(L,Ts[0])
        if Wolff:
            system.mcmc = system.Wolff
        for T in Ts:
            system.reset(T)
            for t in range(0,Nit,Nskip):
                system.mcmc(Nskip)
                plt.cla()
                plt.imshow(system.s, vmin=-1, vmax=+1)
                plt.title(f"T = {T}  Time = {t}")
                plt.pause(.0000001)
        return system

    @staticmethod
    def demo(L : int, T : float, Nit = 1000000, Nskip = 10, Wolff = False):
        system = Ising(L,T)

        if Wolff:
            system.mcmc = system.Wolff

        system.reset(T)
        
        # Set up interactive window:
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.05, bottom=0.05)
        ax.axis('square')

        # Wolff toggle button
        axWolffBut = plt.axes([0.83, 0.05, 0.15, 0.1])
        WolffButton = Button(axWolffBut, 'Metropolis', color='lightgoldenrodyellow', hovercolor='0.975')

        # Temperarure slider
        axT = plt.axes([0.9, 0.25, 0.02, 0.6], facecolor='lightgoldenrodyellow')
        T_slider = Slider(axT, 'T', 0.001, 6.0, orientation='vertical', valinit=T, )

        def T_changed(T):
            system.reset(T)
            fig.canvas.draw_idle()

        T_slider.on_changed(T_changed)

        def toggle_Wolff(mouse_event):
            if WolffButton.label.get_text() == 'Wolff':
                WolffButton.label.set_text('Metropolis')
                system.mcmc = system.metropolis
            elif WolffButton.label.get_text() == 'Metropolis':
                WolffButton.label.set_text('Wolff')
                system.mcmc = system.Wolff
            fig.canvas.draw_idle()

        WolffButton.on_clicked(toggle_Wolff)
    
        for t in range(0,Nit,Nskip):
            system.mcmc(Nskip)
            ax.cla()
            ax.imshow(system.s, vmin=-1, vmax=+1)
            # plt.title(f"T = {T}  Time = {t}")
            plt.pause(.0000001)
        return system

    @staticmethod
    def run_simulation(L : int, Ts : np.iterable, nwarmup : int, nsample : int, filename : str, Wolff = False):
        system = Ising(L,Ts[0])
        if Wolff:
            system.mcmc = system.Wolff
        m = []
        absm = []
        e = []
        chi = []
        cv = []
        binder = []

        with open(filename,'w') as file:
            for T in Ts:
                print(f"{T=}")
                system.reset(T)
                system.run(nwarmup,nsample,1)
                m.append(system.m)
                absm.append(system.absm)
                e.append(system.e)
                chi.append(system.chi)
                cv.append(system.cv)
                binder.append(system.binder)
                system.write(file)

        plt.figure()
        plt.plot(Ts, m)
        plt.plot(Ts, absm)

        plt.figure()
        plt.plot(Ts, e)
        plt.plot(Ts, cv)

        plt.figure()
        plt.plot(Ts, chi)
        plt.figure()
        plt.plot(Ts, binder)

        plt.show()
        return system


# np.random.seed=12345678   # Uncomment for reproducible runs using same seed.
if __name__ == "__main__":

    import sys
    sys.setrecursionlimit(5000) # for Wolff

    L = 32
    nsample = 100000

    Tc = 2.269185314213022

    Ising.demo(L, Tc, nsample, 1, Wolff=False)

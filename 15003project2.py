import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sys


class DarkMatter():

    def __init__(self, N, sigmav):

        # Input parameters
        self.N = N                      # Number of steps
        self.sigmav = sigmav            # Mean Cross section

        # Constants
        self.mxi = 1000.                # GeV, mass
        self.G = 6.674e-11              # Nm^2/kg^2,  Gravitational constant
        self.g_stars = 11/2.               # Degrees of freedom in entropy density
        self.gstar = 28. + 7/8.*90.     # Degrees of freedom for all particle types
        self.g = 2.                     # Degrees of freedom of DM
        self.eps = 0.05                 # +/- difference from DM abundance

        # Arrays
        self.x = np.logspace(0, 3, self.N)
        self.W = np.zeros(self.N)
        self.Y = np.zeros(self.N)
        self.Xsec = np.linspace(1e-14, 1e-7, self.N)


    def yeq(self,x):
        """
        The quantity y = lambda*Y at equilibrium
        """

        mass = self.mxi/(1000.)          # dimensionless
        a = self.sigmav/(1e-10)          # <sigma*v>/1e10GeV^-2
        factor = 9.35e+9*self.g/2.

        yeq = factor*mass*a*np.sqrt(100./self.gstar)*x**(3./2.)*np.exp(-x)
        return yeq


    def Lambda(self):
        """
        Compute Lambda
        """

        factor = 2*np.pi*np.sqrt(90)/(45*np.sqrt(8*np.pi*self.G))
        l = factor*self.g_stars*self.sigmav/np.sqrt(self.gstar)
        return l


    def Weq(self,x):
        """
        The parameter W at equilibrium
        """

        yeq = self.yeq(x)
        Weq = np.log(yeq)
        return Weq



    def dWdx(self, x, W):
        """
        The differential equation to be integrated.
        """

        expfac = np.exp((2*self.Weq(x)-W)) - np.exp(W)
        dWdx = expfac/(x**2)
        return dWdx


    def Integrate(self):
        """
        Integrate the differential equation for the Boltzmann function of Dark
        Matter. Use Forward Euler in the integration.
        """

        x = self.x
        W = self.W

        Weq = self.Weq(x)
        yeq = self.yeq(x)
        W[0] = Weq[0]

        # Integrate
        results = solve_ivp(self.dWdx, (x[0], x[-1]),[W[0]], method='Radau', t_eval=x)
        self.xx = results.t
        W = np.squeeze(results.y)
        for i in range(len(W)):
            # Find the x values when decoupling
            if abs(W[i] - Weq[i]) >= 0.05:
                print i, self.xx[i], W[i]
                break

        self.y = np.exp(W)
        return self.y


    def Ex6(self, Xsection):
        """
        From the given inital cross section values, find out who gives the closest
        Omega_i0 to the measured one from the LambdaCDM model.
        """

        y = self.y
        self.LCDM_Omega = 0.12
        xflist = []
        print '-----------'
        for i in range(len(y)):
            if y[i] <= 0.1*y[0]:
                xflist.append(self.xx[i])

        xf = xflist[0]
        print 'xf =', xf

        self.factor = 1.69*(xf/20.)*np.sqrt(100./self.gstar)
        print 'For LCDM: <sigma v> = ', self.factor*1e-10/0.12
        OmegaDM0 = self.factor*(1e-10/(Xsection))

        print 'Omega_DM0 for different cross sections:'
        for k in range(len(Xsection)):
            print 'Input <sigma v> = %0.1e'%(Xsection[k]),'gives Omega_DM0 = %.3g'%(OmegaDM0[k])

        self.xf = xf
        return OmegaDM0


    def Xsection(self):
        """
        Find the best fitting value of the cross section in the range 1e-14 to
        1e-7, to find the value giving the closest value to dark matter abundance,
        OmegaDM0 = 0.12 +/- 0.05.
        """

        OmegaDM0 = self.factor*(1e-10/self.Xsec)
        BestFit = []; ind = []; Xsec_val = []
        for i in range(self.N):
            if OmegaDM0[i] >= self.LCDM_Omega - self.eps:
                if OmegaDM0[i] <= self.LCDM_Omega + self.eps:
                    BestFit.append(float(OmegaDM0[i]))
                    ind.append(i)
                    Xsec_val.append(float(self.Xsec[i]))

        print '-------------'
        print 'Len of the list with best fitting Omega_DM0:', len(BestFit)

        if len(BestFit) == 1:
            print 'The best value of Dark Matter abundance is: OmegaDM0 =',(BestFit[0])
            print 'The best cross section for this Omega is: <sv> =',(Xsec_val[0])

        elif len(BestFit) == 0:
            print 'No values close to observed Dark Matter abundance,', self.eps

        else:

            under = []; over = []; XSover = []; XSunder = []  # Lists for finding best values
            for j in range(len(BestFit)):
                if BestFit[j] < self.LCDM_Omega:
                    under.append(BestFit[j])
                    XSunder.append(Xsec_val[j])
                if BestFit[j] >= self.LCDM_Omega:
                    over.append(BestFit[j])
                    XSover.append(Xsec_val[j])

            if abs(self.LCDM_Omega - max(under)) > abs(self.LCDM_Omega - min(over)):
                print 'The best value of Dark Matter abundance is Omega_DM0 =', min(over)
                print 'The best cross section for this Omega is: <sv> =',(max(XSover))

            else:
                print 'The best value of Dark Matter abundance is Omega_DM0 =', max(under)
                print 'The best cross section for this Omega is: <sv> =',(min(XSunder))

        print 'Span of values: Omega_DM0 = %g to %g'%(min(over), max(under))
        print 'Span of cross section: %.2e to %.2e GeV^2'%(min(XSover), max(XSunder))

        #### Plot the dark mater abundance ####
        plt.figure('Dark Matter abundance')
        plt.semilogx(self.Xsec, abs(OmegaDM0-self.LCDM_Omega), '-b', label=r'$|\Omega_{DM,0}-0.12|$')
        plt.axhline(y=0.05, color='r', label='margin error = 0.05')

        plt.legend(bbox_to_anchor=(0.01, 0.99), loc=2)
        plt.xlabel(r'$Cross section$')
        plt.ylabel(r'$\Omega_{DM,0}$')
        plt.ylim(-0.05,0.3)
        plt.savefig('DM_abundance.png')

        return BestFit


    def PlotDM(self):
        """
        Plot of the Boltzmann function
        """

        y = self.y
        
        plt.figure('Boltzmann of %.2e'%(self.sigmav))

        plt.loglog(self.xx, y, '-b', label='y(x)')
        plt.loglog(self.xx, self.yeq(self.xx), '-r', label=r'$y_{eq}$')

        plt.legend(bbox_to_anchor=(0.8, 0.99), loc=2)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.ylim(1e-6,1e+11)
        #plt.title(r'Boltzmann function when <$\sigma v$> = %.2e'%(self.sigmav))
        plt.savefig('BoltzmannXS%s.png'%(self.sigmav))




####### Call on the rutine ########

N = int(2e+4)
sigmav = np.array([1e-9, 1e-10, 1e-11])           # GeV^-2


for i in range(len(sigmav)):
    DM = DarkMatter(N, sigmav[i])
    #l = DM.Lambda()
    integrate = DM.Integrate()
    p1 = DM.PlotDM()

ex6 = DM.Ex6(sigmav)
Xsec = DM.Xsection()

plt.show()
print '------------'

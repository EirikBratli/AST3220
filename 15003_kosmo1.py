import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import scipy.fftpack
from sklearn.preprocessing import Imputer

imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
print 'Choose at the bottom of the .py-file what functions to run. default is all functions'


class Project1():

    def __init__(self, N):

        self.n = int(N)
        self.alpha = 1.
        self.iota = 3./2.
        self.kappa = 8*np.pi*6.674e-11 # Nm^2 kg^-2
        self.zi = 2e+7          # redshift at big bang
        h = 0.68                # Measured value of today
        self.H0 = 100*h         # km s^-1 Mpc^-1

        # Convert factors
        self.kmperMpc = 1000/(3.085678e+16*1e+6)
        self.stoyear = (60**2*24.25*365.25*1e+9)
        self.N, self.dN = np.linspace(np.log(1./(1+self.zi)), 0, self.n+1, retstep=True)
        self.z = np.exp(-self.N)-1             # redshift grid
        self.V0 = 1


    def PowPot(self, phi):
        """
        function for the inverse power law potential
        """
        M = 1
        V = M**(4+self.alpha)*phi**(-self.alpha)
        return V

    def ExpPot(self, phi):
        """
        Function for the exponential potential
        """
        V = self.V0*np.exp(-self.kappa*self.iota*phi)
        return V

    def dx1dN(self, x1, x2, x3, l):
        """
        Function for the derviative of x1
        """

        fac = 3*(1+x1**2-x2**2+x3**2/3.)
        return -3*x1 + 0.5*np.sqrt(6)*l*x2**2+0.5*x1*fac

    def dx2dN(self, x1, x2, x3, l):
        """
        Function for the derviative of x2
        """

        fac = 3*(1+x1**2-x2**2+x3**2/3.)
        return -0.5*np.sqrt(6)*l*x1*x2 + 0.5*x2*fac

    def dx3dN(self, x1, x2, x3, l):
        """
        Function for the derviative of x3
        """

        fac = 3*(1+x1**2-x2**2+x3**2/3.)
        return -2*x3 + 0.5*x3*fac

    def dldN(self, l, G, x1):
        """
        Function for the derviative of lambda
        """

        dldN = -np.sqrt(6)*l**2*(G-1)*x1
        return dldN

    def Gpow(self):
        """
        The Gamma factor as a function for the inverse power law potential
        """

        Gamma = 2.0
        return Gamma

    def Gexp(self):
        """
        The Gamma factor as a function for the inverse power law potential
        """

        Gamma = 1.0
        return Gamma


    def Integrate(self, x10, x20, x30, l0,G):
        """
        The integration function, using forward Euler with small steplength,
        return an array with the lambda, x1, x2 and x3 values.
        """
        n = len(self.N)
        X = np.zeros((n, 4))

        X[0,0] = l0
        X[0,1] = x10
        X[0,2] = x20
        X[0,3] = x30

        # Euler:
        for i in range(n-1):

            X[i+1,0] = X[i,0] + self.dN*self.dldN(X[i,0],G,X[i,1])
            X[i+1,1] = X[i,1] + self.dN*self.dx1dN(X[i,1],X[i,2],X[i,3],X[i,0])
            X[i+1,2] = X[i,2] + self.dN*self.dx2dN(X[i,1],X[i,2],X[i,3],X[i,0])
            X[i+1,3] = X[i,3] + self.dN*self.dx3dN(X[i,1],X[i,2],X[i,3],X[i,0])

            if i%500000 == 0:
                print 'i=',i,'-', X[i,:]

        return X, self.z


    def EquationOfState(self,X):
        """
        Compute the equation of state for the given potential
        """

        print 'Calculating equation of state'
        X1 = X[:,1]
        X2 = X[:,2]

        rho_phi = X1**2 + X2**2
        P_phi = X1**2 - X2**2
        omega_phi = P_phi/rho_phi

        return omega_phi


    def Hubble_z(self, X, w_phi):
        """
        Calculate the Hubble parameter as a function of z.
        """

        print 'Calculate Hubble parameter as function of z'

        one = np.ones(self.n+1)            # Array of ones
        zint = np.zeros(self.n+1)           # Arrauy for integration

        zint[0] = 3*(1+w_phi[-1])/(1+self.z[-1])
        self.H02 = self.H0**2

        Omega_r = X[:,3]**2; Omega_phi = X[:,1]**2 + X[:,2]**2
        Omega_m = one - (X[:,1]**2+X[:,2]**2+X[:,3]**2)
        Omega_r0, Omega_m0, Omega_phi0 = Omega_r[-1], Omega_m[-1], Omega_phi[-1]
        self.Omega_pow = [Omega_r, Omega_m, Omega_phi]

        self.dz = -self.dN*np.exp(-self.N)
        for i in range(self.n-1):
            zint[i+1] = zint[i] + self.dz[i]*3*(1+w_phi[i])/(1+self.z[i])

        zterm1 = Omega_r0*(one + self.z)**4
        zterm2 = Omega_m0*(one + self.z)**3
        zterm3 = Omega_phi0*np.exp(zint)
        self.zterms = [zterm1,zterm2,zterm3]
        Hz2 = self.H02*(zterm1 + zterm2 + zterm3)
        print 'H0 = ', np.sqrt(Hz2[-1])
        H = np.sqrt(Hz2)/self.H0; var = self.z
        print 'Normalized Hubble parameter today is: H/H0 = %g'%(H[-1])
        return H, var


    def Hubble_N(self,X,w_phi):
        """
        Calculate the Hubble parameter as a function of N.
        """

        print 'Calculate Hubble parameter as function of N'

        N = np.zeros(self.n+1)
        one = np.ones(self.n+1)
        phiexp = np.zeros(self.n+1)
        phiexp[0] = 3*(1+w_phi[-1])

        Omega_r = X[:,3]**2; Omega_phi = X[:,1]**2 + X[:,2]**2
        Omega_m = one - (X[:,1]**2+X[:,2]**2+X[:,3]**2)
        Omega_r0, Omega_m0, Omega_phi0 = Omega_r[-1], Omega_m[-1], Omega_phi[-1]
        self.Omega_pow = [Omega_r, Omega_m, Omega_phi]

        H02 = self.H0**2
        dz = self.zi/float(self.n+1)
        for i in range(self.n-1):
            dN = -dz/(1+self.z[i])
            N[i] = np.log(1 + self.z[i])
            phiexp[i+1] = phiexp[i] + 3*(1+w_phi[i])*dN

        term1 = Omega_r0*np.exp(-4*N)
        term2 = Omega_m0*np.exp(-3*N)
        term3 = Omega_phi0*np.exp(phiexp)
        self.Nterms = [term1, term2, term3]
        H2 = H02*(term1 + term2 + term3)
        print 'H0 = ', np.sqrt(H2[-1])
        H = np.sqrt(H2/H02); var = N

        print 'Normalized Hubble parameter today is: H/H0 = %g'%(H[-1])
        return H, var


    def Age(self, H, name):
        """
        Calculate the age of the Universe with z as the variable.
        """

        print 'Calcualte the age of the Universe for the %s'%(name)

        dz = self.dN*np.exp(-self.N)
        t0H0z= sum(dz/(H*(1+self.z)))
        t0H0N = sum(self.dN/H)

        dimfactor = 1/(self.kmperMpc)
        print 'For z as variable:'
        print 't0H0=', t0H0z
        print 't0=',(t0H0z/self.H0)/self.kmperMpc/self.stoyear

        print 'For N as variable:'
        print 't0H0=', t0H0N
        print 't0=',(t0H0N/self.H0)/self.kmperMpc/self.stoyear

        print '---------------'
        return t0H0


    def Age_LambdaCDM(self):

        print 'Age of the universe with the LambdaCDM model:'
        Omega_m0 = 0.3; Omega_L0 = 0.7
        H0t0 = 2.*(np.arcsinh(np.sqrt(Omega_L0/Omega_m0)))/(3*np.sqrt(Omega_L0))

        dimfactor = 1/(self.kmperMpc)
        t0 = (H0t0/self.H0)/self.kmperMpc/self.stoyear

        print 'H0t0 = %g'%H0t0              # 0.964099
        print 't0 = %g Gyr'%t0              # 13.7202 Gyr
        return H0t0


    def LuminosityDistance(self, H, eps, name):
        """
        The Hubble parameter is already H/H0. Find the index who gives the value
        closest to 0.83. Then integrate up to that index.
        """

        idx = []
        for ind, val in enumerate(self.z):
            if val <= 0.83:
                idx.append(ind)

        index = min(idx)
        dLz = (1+self.z[index])*sum(self.dz[:index]/H[:index])
        dLN = np.exp(-self.N[index])*sum(np.exp(-self.N[index:])*self.dN/H[index:])

        print 'dL(N) = %g c/H0'%dLN
        print 'dL(z) = %g c/H0'%dLN
        if dLz == dLN:
            print 'Good calculations'
            dL = dLN
        else:
            dL = abs(dLN-dLz)/2.

        print '...........'
        return dL



def plot_Pow():
    name = [r'$\Omega_{m}(z)$',r'$\Omega_{r}(z)$',r'$\Omega_{\phi}(z)$']

    Omega_phi = Xpow[:,1]**2 + Xpow[:,2]**2
    Omega_r = Xpow[:,3]**2
    Omega_m = 1 - Omega_r - Omega_phi
    Omega_i = [Omega_m, Omega_r, Omega_phi]
    for j in range(len(name)):
        plt.figure('IPL')
        plt.semilogx(z, Omega_i[j], label='%s'%name[j])

        plt.xlabel(r'z')
        plt.ylabel(r'Fractional Energy density, $\Omega(z)$')
        plt.title('Inverse Power law Potential')
        plt.legend(bbox_to_anchor=(0.81, 0.65), loc=2, borderaxespad=0.)
        plt.xlim(0,z[0])
        plt.grid(True)
        plt.gca().invert_xaxis()
        plt.savefig('pow_x.png')
    plt.show()

def plot_exp():
    name = [r'$\Omega_{m}(z)$',r'$\Omega_{r}(z)$',r'$\Omega_{\phi}(z)$']

    Omega_phi = Xexp[:,1]**2 + Xexp[:,2]**2
    Omega_r = Xexp[:,3]**2
    Omega_m = 1 - Omega_r - Omega_phi
    Omega_i = [Omega_m, Omega_r, Omega_phi]
    for j in range(len(name)):
        plt.figure('Exp')
        plt.semilogx(z, Omega_i[j], label='%s'%name[j])

        plt.xlabel(r'z')
        plt.ylabel(r'Fractional Energy density, $\Omega(z)$')
        plt.title('Exponential Potential')
        plt.legend(bbox_to_anchor=(0.81, 0.65), loc=2, borderaxespad=0.)
        plt.xlim(0,z[0])
        plt.grid(True)
        plt.gca().invert_xaxis()
        plt.savefig('exp_x.png')

    plt.show()

def plot_eos(omega_phi, name):

    plt.figure('%s'%name)
    plt.semilogx(z,omega_phi,'-b', label=r'$\omega_{\phi}$')

    plt.title(r'Equation of state, %s potential'%(name))
    plt.legend(bbox_to_anchor=(0.85, 0.1), loc=2, borderaxespad=0.)
    plt.xlabel('log(z)')
    plt.ylabel(r'$\omega_{\phi}$')
    plt.gca().invert_xaxis()
    plt.savefig('EoS_%s.png'%(name))
    plt.show()

def plot_Hubble(H, var, name):

    plt.figure('Hubble')
    plt.plot(var, H, label=r'$H(z)/H_0$')
    plt.xlabel('z')
    plt.ylabel(r'$H/H_0$')
    plt.legend(bbox_to_anchor=(0.78, 0.5), loc=2, borderaxespad=0.)
    plt.title('Hubble parameter for %s potential'%name)
    plt.gca().invert_xaxis()
    plt.savefig('Hubbleparam_Pot_%s.png'%(name))
    plt.show()


########################################
project1 = Project1(5e+6)

# Inverse Power Law:
Gpow = project1.Gpow()
Xpow,z = project1.Integrate(5e-5,1e-8,0.9999,1e+9, Gpow)
w_phi_pow = project1.EquationOfState(Xpow)
Hpow,varpow = project1.Hubble_z(Xpow,w_phi_pow)
Age_pow = project1.Age(Hpow, 'IPL potential')
dL = project1.LuminosityDistance(Hpow, 1e-5,'inverse power')

p_pow = plot_Pow()
p_eos = plot_eos(w_phi_pow, 'IPL')
p_Hpow = plot_Hubble(Hpow,varpow, 'IPL')


# Exponential Potential:
Gexp = project1.Gexp()
Xexp,z = project1.Integrate(0, 5e-13,0.9999, 1.5, Gexp)
w_phi_exp = project1.EquationOfState(Xexp)
Hexp,varexp = project1.Hubble_z(Xexp,w_phi_exp)
Age_exp = project1.Age(Hexp, 'exp potential')
dL = project1.LuminosityDistance(Hexp, 1e-5,'exponential')

p_exp = plot_exp()
p_eos = plot_eos(w_phi_exp, 'Exp')
p_Hexp = plot_Hubble(Hexp,varexp, 'Exp')

Age_LCDM = project1.Age_LambdaCDM()

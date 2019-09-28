"""Classes and functions for accessing and manipulating tabulated EOS data."""
### Module for accessing and manipulating tabulated EOS data
### STS 09/2019
###
##
import numpy as np
#
#     ANEOS KPA FLAG
#                                TABLE          ANEOS
#     KPAQQ=STATE INDICATOR      =1, 1p    =1, 1p    (eos without melt)
#                                =2, 2p lv =2, 2p liquid/solid plus vapor
#                                          =4, 1p solid  (eos with melt)
#                                          =5, 2p melt   (eos with melt)
#                                          =6, 1p liquid (eos with melt)
#                                =-1 bad value of temperature
#                                =-2 bad value of density
#                                =-3 bad value of material number
#
class EOShugoniot:
    """Class for Hugoniot array from extEOStable."""	
    def __init__(self):
        self.NH = 0
        self.rho = np.zeros(self.NH)   
        self.T = np.zeros(self.NH)   
        self.P = np.zeros(self.NH)   
        self.U = np.zeros(self.NH)   
        self.S = np.zeros(self.NH)   
        self.up = np.zeros(self.NH)   
        self.us = np.zeros(self.NH)
        self.units = ''
#
class EOSvaporcurve:
    """Class for vapor curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rv = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Pv = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Uv = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Sv = np.zeros(self.NT)
        self.units = ''
#
class EOSmeltcurve:
    """Class for melt curve from ANEOS."""	
    def __init__(self):
        self.NT = 0
        self.NV = 0
        self.T = np.zeros(self.NT)  
        self.rl = np.zeros(self.NT)  
        self.rs = np.zeros(self.NT)  
        self.Pl = np.zeros(self.NT)  
        self.Ps = np.zeros(self.NT)  
        self.Ul = np.zeros(self.NT)  
        self.Us = np.zeros(self.NT)  
        self.Sl = np.zeros(self.NT)  
        self.Ss = np.zeros(self.NT)
        self.units = ''
#
class EOS1barcurve:
    """Class for 1bar curve from the EOS."""	
    def __init__(self):
        self.NT  = 0
        self.S   = np.zeros(self.NT)  
        self.T   = np.zeros(self.NT)  
        self.rho = np.zeros(self.NT)
        self.units = ''
#
class EOScriticalpoint:
    """Class for critical point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.S   = 0  
        self.T   = 0 
        self.rho = 0
        self.U   = 0
        self.units = ''
#
class EOStriplepoint:
    """Class for triple point state from the EOS."""	
    def __init__(self):
        self.P   = 0
        self.T   = 0 
        self.units = ''
#
class EOSaneoshugoniot:
    """Class for Hugoniot calculated in ANEOS."""	
    def __init__(self):
        self.ND  = 0
        self.NV  = 0
        #self.all = np.zeros((self.ND,self.NV))
        self.rho = 0
        self.T   = 0
        self.P   = 0
        self.U   = 0
        self.S   = 0
        self.us  = 0
        self.up  = 0
        self.units = ''
#
class extEOStable:
    """Class for accessing EXTENDED SESAME-STYLE EOS tables output from ANEOS"""
    def __init__(self):
        self.ND  = 0 # integer; number of density points in grid
        self.NT  = 0 # integer; number of temperature points in grid
        self.rho = np.zeros(self.ND)          # g/cm3, density values
        self.T   = np.zeros(self.NT)          # K, temperature values
        self.P   = np.zeros(self.ND*self.NT)  # GPA, pressure(T,rho)
        self.U   = np.zeros(self.ND*self.NT)  # MJ/kg, sp. internal energy(T,rho)
        self.A   = np.zeros(self.ND*self.NT)  # MJ/kg, Helmholtz free energy(T,rho)
        self.S   = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. entropy(T,rho)
        self.cs  = np.zeros(self.ND*self.NT)  # cm/s, sound speed(T,rho)
        self.cv  = np.zeros(self.ND*self.NT)  # MJ/K/kg, sp. heat capacity(T,rho)
        self.KPA = np.zeros(self.ND*self.NT)  # integer, ANEOS KPA flag(T,rho)
        self.MDQ = np.zeros(self.ND*self.NT)  # integer, Model Development Quality Flag(T,rho)
        self.units = ''
        self.hug = EOShugoniot()
        self.vc  = EOSvaporcurve()
        self.mc  = EOSmeltcurve()
        self.cp  = EOScriticalpoint()
        self.tp  = EOStriplepoint()
        self.onebc = EOS1barcurve()
        self.anhug = EOSaneoshugoniot()
        # these are variables needed for the sesame header
        self.MATID   = 0.
        self.DATE    = 0.
        self.VERSION = 0.
        self.FMN     = 0.
        self.FMW     = 0.
        self.R0REF   = 0.
        self.K0REF   = 0.
        self.T0REF   = 0.
        self.P0REF   = 0.
        # variables needed for the ANEOS gamma function
        self.gamma0 = 0.
        self.theta0 = 0.
        self.C24    = 0.
        self.C60    = 0.
        self.C61    = 0.
        self.beta   = 0.
        # model name/version string
        self.MODELNAME = ''

    def loadsesame(self, fname, unitstxt=None):
        """Function for loading EXTENDED SESAME-STYLE EOS table output from ANEOS"""
        data = ([])
        if unitstxt is None:
            self.units = 'Units: rho g/cm3, T K, P GPa, U MJ/kg, A MJ/kg, S MJ/K/kg, cs cm/s, cv MJ/K/kg, KPA flag. 2D arrays are (NT,ND).'
        else:
            self.units = unitstxt
        sesamefile = open(fname,"r")  
        sesamedata=sesamefile.readlines()
        sesamefile.close()
        nskip = 6 # skip standard header to get to the content of the 301 table
        # num.density, num. temps
        tmp = sesamedata[nskip][0:16]
        dlen = float(tmp)
        tmp = sesamedata[nskip][16:32]
        tlen = float(tmp)
        neos = int((dlen*tlen*7.0+dlen+tlen+2.0)/5.0)+1
        #print(dlen,tlen,neos)
        data = np.zeros((neos,5),dtype=float)
        for j in range(nskip,neos+nskip):
            tmp3 = sesamedata[j]
            tmp4 = list(tmp3.split())
            if len(tmp4) < 5:
                lentmp4 = len(tmp4)
                data[j-nskip,0:lentmp4] = np.asarray(tmp4[0:lentmp4])
            else:
                data[j-nskip,:] = np.asarray(tmp4)        
            #print(j,eosarr[j,:])
        #print(data.shape)
        data=np.resize(data,(neos*5))
        #print(data.shape)
        self.ND  = data[0].astype(int)  # now fill the extEOStable class
        self.NT  = data[1].astype(int)
        self.rho = data[2:2+self.ND]
        self.T   = data[2+self.ND : 2+self.ND+self.NT]
        self.P   = data[2+self.ND+self.NT : 2+self.ND+self.NT+self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.U   = data[2+self.ND+self.NT+self.ND*self.NT
                            : 2+self.ND+self.NT+2*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.A   = data[2+self.ND+self.NT+2*self.ND*self.NT
                            : 2+self.ND+self.NT+3*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.S   = data[2+self.ND+self.NT+3*self.ND*self.NT
                            : 2+self.ND+self.NT+4*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cs  = data[2+self.ND+self.NT+4*self.ND*self.NT
                            : 2+self.ND+self.NT+5*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.cv  = data[2+self.ND+self.NT+5*self.ND*self.NT
                            : 2+self.ND+self.NT+6*self.ND*self.NT
                            ].reshape(self.NT,self.ND)
        self.KPA = data[2+self.ND+self.NT+6*self.ND*self.NT
                            : 2+self.ND+self.NT+7*self.ND*self.NT
                            ].reshape(self.NT,self.ND)

    def view(self, q='P', Tlow=None, Thigh=None, rholow=None, rhohigh=None):
        """Function for printing values from EXTENDED SESAME-STYLE EOS table."""
        if Tlow is None:
            Tlow = self.T.min()
        if Thigh is None:
            Thigh = self.T.max()
        if rholow is None:
            rholow = self.rho.min()
        if rhohigh is None:
            rhohigh = self.rho.max()
        print(self.units)
        if q == 'P':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('P:', (self.P[np.logical_and(self.T >= Tlow,
                                                self.T<=Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'U':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('U:', (self.U[np.logical_and(self.T >= Tlow,
                                                self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'A':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('A:', (self.A[np.logical_and(self.T >= Tlow,
                                                self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'S':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho
                                   >= rholow,self.rho<=rhohigh)
                                  ])
            print('S:', (self.S[np.logical_and(self.T >= Tlow,self.T <= Thigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])      
        if q == 'cs':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cs:', (self.cs[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'cv':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cv:', (self.cv[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'KPA':
            print('T:', self.T[np.logical_and(self.T >= Tlow,self.T <= Thigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('KPA:', (self.KPA[np.logical_and(self.T >= Tlow,
                                                  self.T <= Thigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])

    def calchugoniot(self, r0=None, t0=None, pmax=None, writefilename=None):
        """Function for calculating a Hugoniot from EXTENDED SESAME-STYLE EOS table."""
        if r0 is None:
            return 'Must provide r0 and t0.'
        if t0 is None:
            return 'Must provide r0 and t0.'
        if pmax is None:
            pmax=1.E4 # GPa
        self.hug.rho = []
        self.hug.P = []
        self.hug.T = []
        self.hug.U = []
        self.hug.S = []
        self.hug.up = []
        self.hug.us = []
  
        it0 = int(np.round(np.interp(t0,self.T,np.arange(self.NT)))) # uses nearest value if t0 not in array
        ir0 = int(np.round(np.interp(r0,self.rho,np.arange(self.ND)))) # uses nearest value if r0 not in the array
        p0  = self.P[it0,ir0] # GPa
        #print(self.P[it0,ir0])
        e0  = self.U[it0,ir0]#np.interp(p0,self.P[it0,:],self.U[it0,:])
        s0  = self.S[it0,ir0]#np.interp(p0,self.P[it0,:],self.S[it0,:])
        up0 = 0. # no initial particle velocity
        us0 = self.cs[it0,ir0]/1.e5 # cm/s->km/s use sound velocity for initial
        #print(ir0,it0,r0,t0,p0,e0,up0,us0)
        self.hug.rho = np.append(self.hug.rho, self.rho[ir0])
        self.hug.P = np.append(self.hug.P, p0)
        self.hug.T = np.append(self.hug.T, self.T[it0])
        self.hug.U = np.append(self.hug.U, e0)
        self.hug.S = np.append(self.hug.S, s0)
        self.hug.up = np.append(self.hug.up, up0)
        self.hug.us = np.append(self.hug.us, us0)
        
        #for iir in range(ir0+1,self.ND):
        iir=ir0+1
        pnew=p0
        while pnew<pmax:
            ediff =0.5*(self.P[it0::,iir]+p0)*(1./r0-1./self.rho[iir])+e0 -(self.U[it0::,iir])  # MJ/kg
            pnew = np.interp(0.,np.flip(ediff),np.flip(self.P[it0::,iir]))
            tnew = np.interp(0.,np.flip(ediff),np.flip(self.T[it0::]))
            enew = np.interp(0.,np.flip(ediff),np.flip(self.U[it0::,iir]))
            snew = np.interp(0.,np.flip(ediff),np.flip(self.S[it0::,iir]))
            upnew = np.sqrt((pnew-p0)*(1./r0-1./self.rho[iir]))
            usnew = (1./r0)*np.sqrt((pnew-p0)/(1./r0-1./self.rho[iir]))
            #print(self.rho[iir],tnew,pnew,enew,upnew,usnew)
            self.hug.rho = np.append(self.hug.rho, self.rho[iir])
            self.hug.P = np.append(self.hug.P, pnew)
            self.hug.T = np.append(self.hug.T, tnew)
            self.hug.U = np.append(self.hug.U, enew)
            self.hug.S = np.append(self.hug.S, snew)
            self.hug.up = np.append(self.hug.up, upnew)
            self.hug.us = np.append(self.hug.us, usnew)
            iir += 1
        self.hug.NH=len(self.hug.P)
        self.hug.units='units: T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg, Up km/s, Us km/s'

        if writefilename:
            print('Writing Hugoniot to file: ',writefilename)
            hugoniotfile = open(writefilename,"w")  
            hugoniotfile.writelines('  Hugoniot \n') 
            hugoniotfile.writelines('  Temperature,    Density,        Pressure,       Int. Energy,    Sp. Entropy,    Part. Vel.,     Shock Vel. \n') 
            hugoniotfile.writelines('  K,              g/cm3,          GPa,            MJ/kg,          MJ/K/kg,        km/s,           km/s\n') 
            for iih in range(0,self.hug.NH):
                hugoniotfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (
                    self.hug.T[iih],self.hug.rho[iih],self.hug.P[iih],self.hug.U[iih],self.hug.S[iih],self.hug.up[iih],self.hug.us[iih]))
            hugoniotfile.close() 

    def writestdsesame(self, writestdsesfname=None):
        """Write standard Header-201-301 SESAME EOS TABULAR EOS FILE"""
        # write a standard SESAME ascii file
        #     WRITE STANDARD Header-201-301 SESAME FILE
        #     WRITE SESAME 301 TABLE CONTAINS P, E, HFE
        #sesfile = open("NEW-SESAME-STD-NOTENSION.EOSTXT","w")  
        if writestdsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writestdsesfname,"w")  
        #     WRITE SESAME HEADER INFORMATION: EOS matid number, number of words in section
        #     could input matid, date, version with the grid
        # these parameters are set in the cell above that sets up the grid for ANEOS
        # THEY SHOULD MATCH.......
        # These variables are needed for the standard table output
        NWDS=9
        SESNTABLES=2.0
        TABLE1 = 201.0
        TABLE2 = 301.0
        #     5 entries in 201 table
        SESNWDS1=5.0
        #     Number of entries in STANDARD 301 table: 3 variables at each T,rho point
        SESNWDS2=2.+self.ND+self.NT+self.ND*self.NT*3.
        #     HEADER SECTION
        #sesfile.write("%14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e, %14.6e\n" % (antarr[iit],rnew,pnew,enew,snew,upnew,usnew))
        sesfile.write(" INDEX      MATID ={:7d}    NWDS = {:8d}\n".format(int(self.MATID), int(NWDS)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.MATID, self.DATE, self.DATE, self.VERSION, SESNTABLES))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(TABLE1, TABLE2, SESNWDS1, SESNWDS2))
        # 201 SECTION
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE1),int(SESNWDS1)))
        sesfile.write("{:16.8e}{:16.8e}{:16.8e}{:16.8e}{:16.8e}\n".format(self.FMN, self.FMW, self.R0REF, self.K0REF, self.T0REF))
        sesfile.write(" RECORD     TYPE ={:5d}     NWDS = {:8d}\n".format(int(TABLE2),int(SESNWDS2)))
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  pressure array GPa P[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.P[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  specific internal energy array MJ/kg U[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.U[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # Helmholtz free energy array in MJ/kg A[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.A[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the STD SESAME 301 notension table to local directory: ',writestdsesfname)

    def writemdqsesame(self, writemdqsesfname=None):
        """Function to write a sesame 301-style ascii file with the MDQ variable"""
        if writemdqsesfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writemdqsesfname,"w")  
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NT))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     temperature array K
        for j in range(0, int(self.NT)):
            sesfile.write("{:16.8e}".format(self.T[j]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  MDQ Flag[tempindex,dindex]
        for j in range(0,int(self.NT)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.MDQ[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the SESAME TABLE FILE
        sesfile.close() 
        print('Done writing the MDQ Flag as a 301-style table to local directory: ',writemdqsesfname)

    def loadaneos(self, aneosinfname=None, aneosoutfname=None):
        """Function for reading in ANEOS INPUT and OUTPUT FILE DATA into EOS structure."""
        if aneosinfname is None:
            return 'Must provide input file name.'
        if aneosoutfname is None:
            return 'Must provide output file name.'
        # function to gather data from ANEOS input and output files
        # SESAME FILE HEADER INFORMATION MUST BE LOADED INTO THE EOS STRUCTURE BEFORE CALLING THIS FUNCTION
        #
        # READ IN ANEOS INPUT FILE
        aneosinputfile = open(aneosinfname,"r")  
        testin=aneosinputfile.readlines()   # read in the whole ascii file at once because this is fatser
        aneosinputfile.close()
        # gather EOS information from the ANEOS.OUTPUT file
        aneosoutputfile = open(aneosoutfname,"r")  
        testout=aneosoutputfile.readlines() # read everything in at once because this is faster
        aneosoutputfile.close()
        print('Done loading ANEOS files.')

        # THIS CODE PARSES THE ANEOS.OUTPUT FILE INTO ARRAYS FOR USE IN PLOTTING/USING THE EOS
        #print('ANEOS WAS CALLED WITH THE FOLLOWING INPUT, LOADED FROM FILE ',aneosinfname)
        # Gather parameters for the gamma function while printing the ANEOS INPUT FILE
        aneoscount=1
        for i in np.arange(len(testin)):
            if testin[i].find('ANEOS') == 0:
                if aneoscount<9:
                    #print(' '+testin[i-3],testin[i-2],testin[i-1],testin[i])
                    aneoscount=aneoscount+1
                #else:
                    #print(' '+testin[i])
                if testin[i].find('ANEOS2') == 0:
                    tmp=testin[i]
                    nelem=int(tmp[10:20])
                    #print('nelem=',nelem)
                    rho0=float(tmp[30:40])
                    #print('rho0=',rho0)
                    gamma0=float(tmp[70:80])
                    #print('gamma0=',gamma0)
                    theta0=float(tmp[80:90])
                if testin[i].find('ANEOS3') == 0:
                    tmp=testin[i]
                    C24=float(tmp[20:30])/3.
                    #print('C24=',C24)
                if testin[i].find('ANEOS5') == 0:
                    tmp=testin[i]
                    C60=float(tmp[60:70])
                    C61=float(tmp[70:80])
                    #print('C60=',C60)
                if testin[i].find('ANEOS7') == 0:
                    tmp=testin[i]
                    betagamma=float(tmp[70:80])

        # some checks
        #if rho0 != self.R0REF:
        #    print('WARNING: rho0 does not match. STOPPING THIS NOTEBOOK.')
        #    assert(False) # just a way to stop the notebook

        # GUESS A BIG ARRAY SIZE FOR THE PHASE BOUNDARIES AND HUGONIOT IN ANEOS.OUTPUT
        # the melt curve, vapor curve and Hugoniot curves are not fixed length outputs
        nleninit=300
        meltcurve = 0

        print('READING DATA FROM ANEOS OUTPUT FILE ',aneosoutfname)

        # Read in data from the ANEOS.OUTPUT FILE
        imc = -1 # flag for no melt curve in the model
        for i in np.arange(len(testout)):
            if testout[i].find('  Data for ANEOS number') == 0:
                tmp = testout[i+2][0:50]
                eosname = tmp.strip()
            if testout[i] == '  TWO-PHASE BOUNDARIES\n':
                nvc = nleninit
                ivc = i
                vcarrtmp = np.zeros((nvc,12),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    if testout[j+i+4].find(' anphas') == 0:
                        print(testout[j+i+4])
                        vcarrtmp[j,:]=vcarrtmp[j-1,:]
                        j=j+1
                    else:
                        tmp=str.replace(testout[j+i+4],'D','E')
                        tmp3 = tmp[0:157]
                        tmp4 = list(tmp3.split())
                        if (len(tmp4) >0) and (float(tmp4[3]) > 0) and (float(tmp4[4]) > 0): # stop if the pressures become negative on the vapor curve
                            tmp5 = np.asarray(tmp4)
                            vcarrtmp[j,:] = tmp5[:]
                            j=j+1
                        else:
                            flag=1
                vcarr = np.zeros((j,12),dtype=float)
                vcarr[:,:] = vcarrtmp[0:j,:]
            if testout[i] == ' LIQUID/SOLID PHASE CURVE\n':
                nmc = nleninit
                imc = i
                meltcurve=1
                mcarrtmp = np.zeros((nmc,11),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    tmp  = str.replace(testout[j+i+5],'D','E')
                    tmp3 = tmp[0:132]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp5 = np.asarray(tmp4)
                        mcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                mcarr = np.zeros((j,11),dtype=float)
                mcarr[:,:] = mcarrtmp[0:j,:]
            if testout[i] == '   HUGONIOT\n':
                nhc = nleninit
                ihc = i
                hcarrtmp = np.zeros((nhc,9),dtype=float)
                flag=0
                j=0
                while flag == 0:
                    tmp=str.replace(testout[j+i+5],'D','E')
                    tmp3 = tmp[0:109]
                    tmp4 = list(tmp3.split())
                    if len(tmp4) > 0:
                        tmp4[3]='0.0' # this column often gives problems with exponential notation so don't read it
                        tmp5 = np.asarray(tmp4)
                        hcarrtmp[j,:] = tmp5[:]
                        j=j+1
                    else:
                        flag=1
                hcarr = np.zeros((j,9),dtype=float)
                hcarr[:,:] = hcarrtmp[0:j,:]

        # UPDATE THE MAIN EOS STRUCTURE WITH GATHERED INFORMATION
        # Add variables needed to calculate the ANEOS gamma function
        self.gamma0  = gamma0
        self.theta0  = theta0
        self.C24     = C24
        self.C60     = C60
        self.C61     = C61
        self.beta    = betagamma
        #
        # ANEOS.OUTPUT UNITS ARE NOT THE SAME AS THE SESAME TABLE!
        # add the vapor curve to this EOS object extracted from the ANEOS.OUTPUT
        #  TWO-PHASE BOUNDARIES
        #       T         RHOLIQ        RHOVAP        PLIQ         PVAP        ELIQ         EVAP         SLIQ         SVAP        GLIQ         GVAP         PSILIQ      PSIVAP         NTY
        #       K         kg/m**3       kg/m**3       GPa          GPa         J/kg         J/kg        J/kg-K       J/kg-K       J/kg         J/kg
        tmp = vcarr.shape
        #put vapor curve information in nicely named structure
        self.vc.NT  = tmp[0]
        self.vc.T   = vcarr[:,0] # K
        self.vc.rl  = vcarr[:,1]/1.E3 # g/cm3
        self.vc.rv  = vcarr[:,2]/1.E3 # g/cm3
        self.vc.Pl  = vcarr[:,3] # GPa
        self.vc.Pv  = vcarr[:,4] # GPa
        self.vc.Ul  = vcarr[:,5]/1.E6 # MJ/kg
        self.vc.Uv  = vcarr[:,6]/1.E6 # MJ/kg
        self.vc.Sl  = vcarr[:,7]/1.E6 # MJ/K/kg
        self.vc.Sv  = vcarr[:,8]/1.E6 # MJ/K/kg
        self.vc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
        #
        # add the ANEOS Hugoniot to this EOS object extracted from the ANEOS.OUTPUT
        #      RHO          T           P          PC           E           S           V           U       RHO/RHOO  #IT  STATE
        #    kg/m**3        K          GPa        GPa          J/kg      J/kg-K       km/sec      km/sec
        #self.anhug.all = hcarr  # 2D array of Hugoniot variables
        tmp = hcarr.shape
        self.anhug.ND = tmp[0] # number of density points on the Hugoniot
        self.anhug.rho = hcarr[:,0]/1.E3 # g/cm3
        self.anhug.T   = hcarr[:,1] # K
        self.anhug.P   = hcarr[:,2] # GPa
        self.anhug.U   = hcarr[:,4]/1.E6 # MJ/kg
        self.anhug.S   = hcarr[:,5]/1.E6 # MJ/K/kg
        self.anhug.us  = hcarr[:,6] # km/s
        self.anhug.up  = hcarr[:,7] # km/s
        self.anhug.units = 'vars: rho g/cm3, T K, P GPa, U MJ/kg, S MJ/K/kg, Us km/s, Up km/s'
        #
        # Add melt curve to EOS objects if available
        # LIQUID/SOLID PHASE CURVE
        #       T         RLIQ       RSOLID      PLIQ       PSOLID      ELIQ        ESOLID       SLIQ       SOLID        GLIQ       GSOLID        #ITER
        #       K        kg/m**3     kg/m**3      GPa         GPa       J/kg         J/kg        J/kg-K     J/kg-K       J/kg        J/kg
        if meltcurve == 1:
            # put the melt curve information in nicely named structure
            tmp=mcarr.shape
            self.mc.NT  = tmp[0]
            self.mc.T   = mcarr[:,0] # K
            self.mc.rl  = mcarr[:,1]/1.E3 # g/cm3
            self.mc.rs  = mcarr[:,2]/1.E3 # g/cm3
            self.mc.Pl  = mcarr[:,3] # GPa
            self.mc.Ps  = mcarr[:,4] # GPa
            self.mc.Ul  = mcarr[:,5]/1.E6 # MJ/kg
            self.mc.Us  = mcarr[:,6]/1.E6 # MJ/kg
            self.mc.Sl  = mcarr[:,7]/1.E6 # MJ/K/kg
            self.mc.Ss  = mcarr[:,8]/1.E6 # MJ/K/kg
            self.mc.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
            self.tp.T = mcarr[3,0] # K
            self.tp.P = mcarr[3,3] # GPa
            self.tp.units = 'T K, P GPa'
        # put the data for the critical point in the EOS structure for easy access
        self.cp.T   = vcarr[0,0] # K
        self.cp.rho = vcarr[0,1]/1.E3 # g/cm3
        self.cp.P   = vcarr[0,3] # GPa
        self.cp.U   = vcarr[0,5]/1.E3 # MJ/kg 
        self.cp.S   = vcarr[0,7]/1.E3 # MJ/K/kg
        self.cp.units = 'T K, rho g/cm3, P GPa, U MJ/kg, S MJ/K/kg'
        #------------------------------------------------------------------

###########################################################################################
########## GADGET STYLE TABLES
###########################################################################################
class extGADtable:
    """Class for v2019 GADGET EOS tables using SESAME UNITS"""
    def __init__(self):
        """Function to initialize GADGET TABULAR EOS structure"""
        self.ND  = 0 # integer; number of density points in grid
        self.NS  = 0 # integer; number of sp. entropy points in grid
        self.rho = np.zeros(self.ND)          # g/cm3, density values
        self.S   = np.zeros(self.NS)          # MJ/K/kg, sp. entropy values
        self.P   = np.zeros(self.ND*self.NS)  # GPA, pressure(S,rho)
        self.T   = np.zeros(self.ND*self.NS)  # K, temperature(S,rho)
        self.U   = np.zeros(self.ND*self.NS)  # MJ/kg, sp. internal energy(S,rho)
        self.A   = np.zeros(self.ND*self.NS)  # MJ/kg, Helmholtz free energy(S,rho)
        self.cs  = np.zeros(self.ND*self.NS)  # cm/s, sound speed(S,rho)
        self.cv  = np.zeros(self.ND*self.NS)  # MJ/K/kg, sp. heat capacity(S,rho)
        self.KPA = np.zeros(self.ND*self.NS)  # integer, ANEOS KPA flag(S,rho)
        self.MDQ = np.zeros(self.ND*self.NS)  # integer, Model Development Quality Flag(S,rho)
        self.units = 'Units: g/cm3, MJ/K/kg, GPa, K, MJ/kg, MJ/kg, cm/s, MJ/K/kg, KPA flag. 2D arrays are (NS,ND).'
        self.MODELNAME = ''
    #
    def view(self, q='P', Slow=None, Shigh=None, rholow=None, rhohigh=None):
        """Function to print variables from GADGET EOS table"""
        if Slow is None:
            Slow = self.S.min()
        if Shigh is None:
            Shigh = self.S.max()
        if rholow is None:
            rholow = self.rho.min()
        if rhohigh is None:
            rhohigh = self.rho.max()
        if q == 'P':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('P:', (self.P[np.logical_and(self.S >= Slow,
                                                self.S<=Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'U':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('U:', (self.U[np.logical_and(self.S >= Slow,
                                                self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'A':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho <= rhohigh)
                                  ])
            print('A:', (self.A[np.logical_and(self.S >= Slow,
                                                self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])
        if q == 'T':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho
                                   >= rholow,self.rho<=rhohigh)
                                  ])
            print('T:', (self.T[np.logical_and(self.S >= Slow,self.S <= Shigh)
                               ])[:, np.logical_and(self.rho >= rholow,
                                                     self.rho <= rhohigh)
                                 ])      
        if q == 'cs':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cs:', (self.cs[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'cv':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('cv:', (self.cv[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
        if q == 'KPA':
            print('S:', self.S[np.logical_and(self.S >= Slow,self.S <= Shigh)])
            print('rho:', self.rho[np.logical_and(self.rho >= rholow,
                                                   self.rho<=rhohigh)])
            print('KPA:', (self.KPA[np.logical_and(self.S >= Slow,
                                                  self.S <= Shigh)
                                 ])[:, np.logical_and(self.rho >= rholow,
                                                       self.rho <= rhohigh)
                                   ])
    #
    def writeextgadget(self, writeextgadgetfname=None):
        """Function to write an extended gadget2 EOS ascii file in 301 format"""
        #sesfile = open("NEW-GADGET-EXT-NOTENSION.EOSTXT","w")  
        if writeextgadgetfname is None:
            print('Please provide a file name.')
            exit(0)
        sesfile = open(writeextgadgetfname,"w")  
        sesfile.write("{:16.8e}{:16.8e}".format(self.ND, self.NS))
        STYLE=2
        #     density array g/cm3
        for k in range(0, int(self.ND)):
            sesfile.write("{:16.8e}".format(self.rho[k]))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #     sp. entropy array MJ/K/kg -> erg/K/g
        for j in range(0, int(self.NS)):
            sesfile.write("{:16.8e}".format(self.S[j]*1.E10))
            STYLE=STYLE+1
            if (np.mod(STYLE,5) == 0):
                sesfile.write("\n")
        #  pressure array GPa -> dynes/cm2 P[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.P[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  temperature array K T[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.T[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  specific internal energy array MJ/kg -> erg/g U[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.U[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        #  sound speed array cm/s cs[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cs[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # Helmholtz free energy array in MJ/kg -> ergs/g A[tempindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.A[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # sp. heat capacity array MJ/K/kg -> erg/K/g cv[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.cv[j,k]*1.E10))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # KPA array flag k[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.KPA[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # MDQ array flag MDQ[sindex,dindex]
        for j in range(0,int(self.NS)):
            for k in range(0,int(self.ND)):
                sesfile.write("{:16.8e}".format(self.MDQ[j,k]))
                STYLE=STYLE+1
                if (np.mod(STYLE,5) == 0):
                    sesfile.write("\n")
        # close the EXTENDED GADGET TABLE FILE
        sesfile.close() 
        print('Done writing the EXTENDED GADGET TABLE FILE with notension to local directory: ',writeextgadgetfname)
#
### END MODULE ###

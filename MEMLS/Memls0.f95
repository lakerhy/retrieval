!--------------------------------------------------------------------------------
!     Interface to Memls written in Matlab code
!     1)  Variable input
!     2)  Fixed input and execution
!     3)  Output
!     Surface emissivity for an array of frequencies and corresponding zenith angles
!     (1)   The common effective temperature required by Rttov for each instrument is taken at the highest value
!           Rttov = frequency independent = grey body
!     (2)   Emissivity is compensated to leave brightness temperature unaffected
!           Storage is [V-pol freq.1, H-pol freq.1, .. , H-pol freq.N]
!--------------------------------------------------------------------------------
      Subroutine Memls0                   (fil,frq0,Nin0,Nin,Nch0,Nch,Nch1,Npar,scene,ice,y,emi,Teff)
! --- Input and internal
      Character(len=2)                    :: ice                                 ! Type of ice
      Character(len=29),dimension(2)      :: fil
      Character(len=200),dimension(2)     :: string
      Integer,Parameter                   :: Ndat=50,Nlayer=11
      Integer                             :: scene,Nin0,Nin,Nch0,Nch1,Npar,Nlayer0,i,j,k,l,p,cat
      Integer,dimension(Nin0)             :: Nch
      Real(kind=8),parameter              :: h=6.6260693d-34,kB=1.3806505d-23, & ! Planck's and Stefan-Boltzmann's constants [J/Hz,J/K]
     &                                       Tsea=271.35d0,ks=0.31d0,k0=2.034d0, & ! (***), [K,W/m/K]
     &                                       beta=0.117d0,smax=0.1d0             ! [W/m/per mil,K/cm] 
      Real(kind=8),dimension(Nin0*Nch0)   :: frq0
      Real(kind=8),dimension(Npar)        :: y                                   ! Geophysical parameters
! --- Memls interface
      Real(kind=8)                        :: ki,T0(2),slope(3),sum0,val0,val(3)
      Real(kind=8),dimension(Nlayer,Nch0) :: Kabs
      Real(kind=8),dimension(Nlayer)      :: Si,T,d,rho,pci,sal,w,dTau,Tau,dsum  ! Memls input profile and output optical depth
! --- Memls output
      Character(len=3),dimension(4)       :: var=(/'fre','gai','Tbh','Tbv'/)     ! Quantity
      Real(kind=8),dimension(Nin*Nch0)    :: Teff0
      Real(kind=8),dimension(2,Nin*Nch0)  :: Tb0,emi0                            ! First index 1 = H-pol and 2 = V-pol
      Real(kind=8),dimension(Nin0,2,Nch0) :: emi1                                ! Emissivity (1=H-pol, 2=V-pol)
! --- Subroutine output
      Real(kind=8),dimension(Nin0)        :: Teff                                ! Effective temperature lowest frequency mixed pol
      Real(kind=8),dimension(Nin0,Nch0)   :: emi                                 ! Emissivity

! -------------------------------------------------------------------------------
!     1. Variable input
! -------------------------------------------------------------------------------
!     Ice:
!     Thickness and salinity transferred
      d(1)=y(7)
      sal(1)=y(9)
! --- Density: One value taken for FY as nilas and MY as below the waterline (10 % freeboard)
      If (ice == 'FY') then
         rho(1)=920d0
      Else
!! RASMUS: OK?
         rho(1)=920d0*0.9d0+900d0*0.1d0
      End if
! --- Volumetric water content [%]: volume of brine [Ulaby et al., 1986, E71] in Vb.m (also Zubov?)
      T0(1)=y(8)-273.15d0
!If (ice == 'FY') then
!w(1)=0d0
!Else
      If (Dabs(T0(1)) < 1d0) the
         w(1)=1d-3*y(9)*49.717d0
      Else
         w(1)=1d-3*y(9)*(-49.185d0/y(8)+0.532d0)
      End if
!End if
!---  Correlation length for scattering agent
!     FY: Scaling by volume of brine yields the desired variability
!         Ma[umlaut]tzler, Relation between grain-size and correlation length of snow, J.Glaciology 48(162), 461-6, 2002
!     MY: Air bubbles entrapped primarily in the top 25 cm
      If (ice == 'FY') then
         pci(1)=0.8d0*((3d0*w(1))/(4d0*3.14d0))**0.333d0
      Else
         If (d(1) >= 0.25d0) then
            pci(1)=0.35d0
         Else
            pci(1)=(0.35d0*0.25d0+0.25d0*(d(1)-0.25d0))/y(7)
         End if
      End if
! --- Snow:
!     Discretized into numerical layers. Thickness, density and correlation length transferred, no salinity or water content
      Do i=2,Nlayer
         d(i)=y(3)/(Nlayer-1)
         sal(i)=0d0
         w(i)=0d0
         rho(i)=y(5)
         pci(i)=y(6)
      End do

! --- Range of validity: (A) Temperature. Memls
      If (y(1) > 271.15d0) then
         Write (1,*) 'Temperature profile adjusted to Memls [-43.2:-2]C = down',y(1)
         y(1)=271.15d0
      Else if (y(1) < 229.95d0) then
         Write (1,*) 'Temperature profile adjusted to Memls [-43.2:-2]C = up  ',y(1)
         y(1)=229.95d0
      End if

! --- Scene: Ice & snow temperature
!     Degrees of freedom = air-sea thermodynamic equilibrium (0) or piecewise linear bended at interface as continuous (1) or not (2)
      If (scene == 0) then
! ---    [Untersteiner, 1986; Bitz CM and WH Lipscomb, 1999: An energy-conserving thermodynamic model of sea ice, JGR 104, 15669-77] 
         T(1)=Tsea
         Do i=1,Ndat
! ---       Heat conductivity of ice when ignoring frost flowers and, in snow, salinity
!!          Multiple salty layers: Identify ki,min 
            ki=k0+(beta*sal(1))/(T(1)+273.15d0)
! ---       Ratio between the temperature gradients in ice and snow
            slope(3)=ki/ks*d(1)/(sum(d)-d(1))
! ---       Temperature at ice-snow inferface
            T0(1)=(y(1)+slope(3)*Tsea)/(1d0+slope(3))
! ---       The heat gradient is insulated by the snow cover from extreme synoptic weather systems
!           [K/cm] positive upwards
            slope(1)=(T0(1)-Tsea)/d(1)
            If (Dabs(slope(1)) > smax) then
               slope(1)=slope(1)/Dabs(slope(1))*smax
               T0(1)=Tsea+slope(1)*d(1)
            End if
            slope(2)=(y(1)-T0(1))/(sum(d)-d(1))
! ---       The fact that the insulative character of snow might force this portion of the profile away from linearity is ignored
!           If (Dabs(slope(2)) > ki*smax)                  slope(2)=slope(2)/Dabs(slope(2))*ki*smax
! ---       The piecewise linear temperature profile is iterated up to Ndat times
            T(1)=Tsea+slope(1)*d(1)/2d0
            T(2)=T0(1)+slope(2)*d(2)/2d0
            Do j=3,Nlayer
               T(j)=T(j-1)+slope(2)*(d(j)+d(j-1))/2d0
            End do
            If (i > 1 .and. Dabs(T0(1)-T0(2)) < 0.01d0)    Go to 10
            T0(2)=T0(1)
         End do
   10    Continue
      Else
! ---    Piecewise linear to interface with continuity across imposed or not
!        The latter serves climatological averaging of numerous profiles. The harmonic-mean conductivity of two layers is
!        a non-linearity that is letting the varying depths push the average away from continuity.
         T0(1)=2d0*y(4)-y(1)
         T0(2)=2d0*y(8)-Tsea
         If (scene == 1) then
            T0(1)=(T0(1)+T0(2))/2d0
            T0(2)=T0(1)
         End if
         T(1)=(T0(2)+Tsea)/2d0
         slope(1)=(y(1)-T0(1))/y(3)
         Do i=2,Nlayer
            T(i)=T0(1)+slope(1)*(i-1.5d0)
         End do
      End if
! --- Update bulk temperature of snow and ice [K]
!     Linearity of profile implies they can be taken as the average of the ice interface and either air or sea
      y(4)=(y(1)+T0(1))/2d0
      y(8)=(Tsea+T0(2))/2d0

! -------------------------------------------------------------------------------
!     2. Fixed input and execution
! -------------------------------------------------------------------------------
!     Index numerical layers as ice below ten layers of snow
      Si=(/1,0,0,0,0,0,0,0,0,0,0/)
! --- Input fil(1)    : frequency [GHz]
!           fil(2)    : incidence angle [deg]
!           syntes.txt: level indexed upwards
!     type            : ice first year = 3 or multiyear = 4
      If (ice == 'FY') then
         string(1)=       "octave -i -q > Input/Memls/Memls-out.txt << END \n&
    &    icemain('Input/Memls/Freq-Memls-in.txt','Input/Memls/Angl-Memls-in.txt','Syntes.txt',3)"C
      Else
         string(1)=       "octave -i -q > Input/Memls/Memls-out.txt << END \n&
    &    icemain('Input/Memls/Freq-Memls-in.txt','Input/Memls/Angl-Memls-in.txt','Syntes.txt',4)"C
      End if
!! Alternativ (direkte i terminal - ej fra Fortran)
!octave --eval "icemain('Input/Memls/Freq-Memls-in.txt','Input/Memls/Angl-Memls-in.txt','Syntes.txt',0,0,0,271.35,3)" -q > 'Input/Memls/Memls-out.txt'
!string(1)='octave --eval '//'"icemain()"'//' -q > Input/Memls/Memls-out.txt'

! --- Input sequence: layer #, temp [K], vol.water content, density [kg/m3], thickness [cm], cor.length [mm], salinity [per mil], ice/snow[1/0]
      Open                (unit=4,file='Syntes.txt',status='unknown')
      Do i=1,Nlayer
         Write            (4,*)                         i,T(i),w(i),rho(i),d(i),pci(i),sal(i),Si(i)
      End do
      Close               (unit=4)
! --- Call Memls using Icemain as interface
      status=             system(string(1))

! -------------------------------------------------------------------------------
!     3. Output
! -------------------------------------------------------------------------------
!     Read frequency, absorption coefficient each layer, and brightness temperature each polarisation
!     Order of frequencies determined in the mapping above
      Kabs=0d0
      Tb0=0d0
      Open                (unit=4,file='Input/Memls/Memls-out.txt',status='old')
      Do i=1,Nch1
         Do j=1,4
            cat=0
! ---       Expansion to Ndat spawns no print statements in any of Memls' 27 subroutines
            Do k=1,Ndat
               Read       (4,'(a40)',End=13)           string(2)
               cat=       Index(string(2),var(j),back=.true.)
               If (cat > 0) then
                  If      (j == 1) then
! ---                Statement on # incoherent layers
!                    Not necessary to reduce # layers with those that are incoherent (****)
                     Do l=1,Ndat
                        Read (4,'(a40)',End=13)        string(2)
                        cat=Index(string(2),'coherent layers',back=.true.)
                        If (cat > 0) then
!! (****) Nlayer0: Slet l¾sning eller brug?
! ---                      Written 'no coherent layers detected:' or '1 coherent layers of 11 detected:'
                           If (string(2)(1:2) == 'no') then
                              Nlayer0=Nlayer
                           Else
                              Read (string(2)(1:cat-1),*) val0
                              Nlayer0=val0
                           End if
                           Go to 13
                        End if
                     End do
                  Else if (j == 2) then
                     Do l=1,Nlayer+1
                        If (l > 1) Read (4,*,End=13)   Kabs(l-1,i)
                     End do
                  Else
                     Read (string(2)(cat+6:40),*)      Tb0(j-2,i)
                  End if
                  Go to 13
               End if
            End do
   13       Continue
         End do
      End do
      Close               (unit=4)

! --- Calculations
      Teff0=0d0
      emi0=0d0
      Do i=1,Nch1
! ---    Incremental optical depth [m^-1]
         Do j=Nlayer,1,-1
            dTau(j)=Kabs(j,i)*d(j)/100d0
         End do
! ---    Accumulate layers above. Indexing progresses upward
         Tau=0d0
         Do j=Nlayer,1,-1
            If (j < Nlayer) then
               Do k=Nlayer,j+1,-1
                  Tau(j)=Tau(j)+dTau(k)
               End do
            End if
         End do
         sum0=0d0
         Do j=Nlayer,1,-1
            dsum(j)=Dexp(-Tau(j))*dTau(j)
            sum0=sum0+dsum(j)
         End do
! ---    Effective temperature as gauged by the sensor, the equivalent of skin temperature
         Do j=Nlayer,1,-1
            Teff0(i)=Teff0(i)+T(j)*dsum(j)
         End do
         Teff0(i)=Teff0(i)/sum0
! ---    Emissivity (*)
         Do j=1,3
            If (j == 3) then
               val0=Teff0(i)
            Else
               val0=Tb0(j,i)
            End if
            val(j)=Dexp((h*frq0(i)*1d9)/(kB*val0))-1d0
         End do
         Do j=1,2
            emi0(j,i)=val(3)/val(j)
         End do
      End do

! --- Transformation to Rttov
      Teff=0d0
      emi1=0d0
      emi=0d0
      k=1
      Do i=1,Nin
! ---    Effective temperature has to be held constant for the given instrument
!        Selecting the largest value ensures all emissivities remain below unity when compensating below
         Teff(i)=Teff0(k)
         Do j=1,Nch(i)
            l=k+j-1
            If (Teff0(l) > Teff(i))                    Teff(i)=Teff0(l)
         End do
! ---    Compensate emissivity at all frequencies
!        Linearisation with the Rayleigh-Jeans approximation deviates less than a per cent from (*)
         Do j=1,Nch(i)
            l=k+j-1
            Do p=1,2
               emi1(i,p,j)=emi0(p,l)*Teff0(l)/Teff(i)
            End do
         End do
         k=k+Nch(i)
! ---    Merge the last two indexes and reverse to V- before H-pol
         Do j=1,Nch(i)
            emi(i,2*j-1)=emi1(i,2,j)
            emi(i,2*j)  =emi1(i,1,j)
         End do
! ---    Range of validity: (B) Emissivity. Got to be physical
         Do j=1,2*Nch(i)
            If (emi(i,j) < 0d0) then
Write (1,*) 'Emissivity adjusted up  ',i,j,emi(i,j)
               emi(i,j)=0d0
            Else if (emi(i,j) > 1d0) then
Write (1,*) 'Emissivity adjusted down',i,j,emi(i,j)
               emi(i,j)=1d0
            End if
         End do
      End do
      End subroutine

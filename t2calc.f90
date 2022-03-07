program t2calc
implicit none

integer           :: i,j,k
double precision  :: tauR, tauR1, tauR2, tauT, tauT1, tauT2, g0R, g0T
double precision  :: convG0, convtau
double precision  :: invT2R,T2R,invT2T,T2T,invT2,T2,T2tmp
character(LEN=50) :: inputfile(3),outputfile

call getarg(1,inputfile(1))
call getarg(2,inputfile(2))
call getarg(3,inputfile(3))
call getarg(4,outputfile)

convG0=1.06679D11
convtau=1E-15

open(10,file=inputfile(1)) !INPUT G0s
read(10,*) g0R,g0T
print*, 'g0R,g0T ',g0R,g0T
close(10)

open(10,file=inputfile(2)) !INPUT tauR
read(10,*) tauR1,tauR2 
print*,'tauR1,tauR2 ',tauR1,tauR2
close(10)

open(10,file=inputfile(3)) !INPUT tauT
read(10,*) tauT1,tauT2
print*,'tauT1,tauT2 ',tauT1,tauT2
close(10)

if(tauR1>tauR2)then
 tauR=tauR1
else
 tauR=tauR2
endif

if(tauT1>tauT2)then
 tauT=tauT1
else
 tauT=tauT2
endif

print*,'tauR,tauT ',tauR,tauT

g0R=g0R*convG0
g0T=g0T*convG0
tauR=tauR*convtau
tauT=tauT*convtau

invT2R=10*g0R*tauR
invT2T=10*g0T*tauT
invT2=invT2R+invT2T
T2=1/invT2
T2R=1/invT2R
T2T=1/invT2T
T2tmp=(T2R*T2T)/(T2R+T2T)

write(*,'(A20,4F12.5)')'invT2,invT2R,invT2T ',invT2,invT2R,invT2T
write(*,'(A12,4F12.5)')'T2,T2tmp,T2R,T2T ',T2,T2tmp,T2R,T2T
open(20,file=outputfile)
write(20,'(3F12.5)')T2,T2R,T2T

end program t2calc


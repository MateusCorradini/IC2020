!*********************************************************************
!		Programa Bloch_sol
!
! Calcula os autovalores na primeira Zona de Brillouin de uma cadeia
! unidimensional monoatômica via método das Somas de Bloch.
!
!*********************************************************************

program bloch_sol
implicit none

	
!------------------ Variáveis a serem utilizadas ---------------------

integer  ::   i, mlim, m, nat
real*8, parameter :: pi = 4.d0*atan(1.0)
real*8   ::   alfa(2), beta(3), a, k
real*8   ::   theta1, theta2
real*8   ::   en, ep

open(9,file="autovalores.dat")

!------------------- Lendo os Parâmetros ------------------------------

print*
write(*,*) ' Entre com o Numero de átomos na cadeia: '
read(*,*) nat
print*
write(*,*) ' Entre com o alpha 1 '
read(*,*) alfa(1)
print*
write(*,*) ' Entre com o alpha 2 '
read(*,*) alfa(2)
print*
write(*,*) ' Entre com o Beta 1 '
read(*,*) beta(1)
print*
write(*,*) ' Entre com o Beta 2 '
read(*,*) beta(2)
print*
write(*,*) ' Entre com o Beta 12 '
read(*,*) beta(3)
print*
write(*,*) ' Entre com o parêmetro de rede a '
read(*,*) a
print*



!------------------ Efetuando os Cálculos --------------------------


mlim = nat/2
m = -mlim

do  i = -mlim, mlim
 	
	k = (2.d0*pi*m)/(nat*a)

	theta1 = alfa(1) + (2.*beta(1)*cos(k*a))
	theta2 = alfa(2) + (2.*beta(2)*cos(k*a)) 
	
	ep = Epos(theta1, theta2, beta(3), k, a)
	en = Eneg(theta1, theta2, beta(3), k, a)
	
	m = m + 1

	write(9,*) k, ep, en

enddo

endfile 9

contains

	real FUNCTION Epos(t1,t2,b12,k,a)  ! Calcula a banda superior
	
	IMPLICIT NONE
	real*8, intent(IN)  ::   t1, t2, b12, k, a
	real*8 :: epo	
	
	epo = 0.5d0*(t1+t2) + 0.5d0*(dsqrt((t1+t2)**2 - 4.0*(t1*t2 - 4.0*(b12**2)*(cos(k*a)**2))))
	
	Epos = epo

	return

	END FUNCTION Epos

       real FUNCTION Eneg(t1,t2,b12,k,a) ! Calcula a banda inferior
	
	IMPLICIT NONE
	real*8, intent(IN)  ::   t1, t2, b12, k, a
	real*8 :: ene

	
	ene = 0.5d0*(t1+t2) - 0.5d0*(dsqrt((t1+t2)**2 - 4.0*(t1*t2 - 4.0*(b12**2)*(cos(k*a)**2))))
 
	Eneg = ene

		return

	END FUNCTION Eneg




end program bloch_sol

!> @file MD.f90 
!>
!> @brief FUNCTIONS: MD 
!>
!> @details A simple NVT MD simulation program.
!> Use Lennard-Jones Potential, 
!> Velocity Verlet Algorithm, 
!> Andersen Thermostat.
!>
!> @author Wenqiang Li 
!>
!>
!> @param Temperture  Temperture, unit of K
!> @param dt          Time interval, unit of s 
!> @param tmax        Max simulation time, unit of s         
!> @param L           Cell length, unit of m
!> @param rc          Cut off, unit of m
!> @param kb          Boltzmann constant, unit of J/K
!> @param m           Mass of Ne, unit of kg
!> @param sigma       Lennard-Jones-sigma
!> @param e           Lennard-Jones-e
!> @param freq        Andersen thermostat collision frequence
!> @param npart       Number of particle
!> @param nsample     Sampling step lengths 
!> @param rc2=rc**2   Sqrt of Cut off, unit of m
 

    
     module constants   


    real(8),parameter::Temperture=275.0    !< Temperture  Temperture, unit of K !
    real(8),parameter::dt=1.0E-16          !< Time interval, unit of s !
    real(8),parameter::tmax=1.0E-12        !< Max simulation time, unit of s 
    real(8),parameter::L=1.0E-9            !< Cell length, unit of m
    real(8),parameter::rc=6.85E-10         !< Cut off, unit of m
    real(8),parameter::kb=1.38E-23         !< Boltzmann constant, unit of J/K
    real(8),parameter::m=3.35E-26          !< Mass of Ne, unit of kg
    real(8),parameter::sigma=2.74E-10      !< Lennard-Jones-sigma
    real(8),parameter::e=5.0E-22           !< Lennard-Jones-e
    real(8),parameter::freq=0.01           !< Andersen thermostat collision frequence
    integer,parameter::npart=100           !< Number of particle          
    integer,parameter::nsample=10          !< Sampling step lengths 
    real(8),parameter::rc2=rc**2           !< Sqrt of Cut off, unit of m
    
    real(8),parameter::pi=3.141592657      !< constant pi

    contains

!> generate inital position and velocity     
    subroutine inital(x,v,a)
    implicit none
    
    integer::i,j,ix,iy,iz,nl,nn
    real(8)::vv1,vv2,vv,vv_sum=0.0
    real(8)::ll 
    
    real(8),intent(out):: x(npart,3),v(npart,3),a(1,3)

    nn=floor(npart**(1.0/3.0))+1
    
    ll=L/real(nn)
    
   i=1
   do iz=1,nn
     do iy=1,nn
        do ix=1,nn
            if (i > npart) exit
            x(i,1)=ix*ll
            x(i,2)=iy*ll
            x(i,3)=iz*ll
            i=i+1
        enddo
     enddo
   enddo
  
    
    call random_seed()
    do j=1,3
        do i=1,npart
            call RANDOM_NUMBER(vv1)
            call RANDOM_NUMBER(vv2)
            vv=sqrt(-2*log(vv1))*cos(2*pi*vv2)
            v(i,j)=vv
            vv_sum=vv_sum+vv**2
        enddo
        a(1,j)=sqrt(npart*kb*Temperture/m/vv_sum)
        v(:,j)=v(:,j)*a(1,j)
    enddo
    
    return

    end subroutine
    
!> calculate force and potential energy of particles
    subroutine cal_F(x,f,E_V)
    implicit none
    
    integer::i,j,k
    real(8)::d(1,3),r,ff(npart,npart,3),f_sum=0.0,fx_t=0.0,fy_t=0.0,fz_t=0.0
    real(8),intent(in)::x(npart,3)
    real(8),intent(out)::f(npart,3),E_V
    
    E_V=0.0
    
    Do i=1,npart-1
        Do j=i+1,npart
            d(1,:)=x(i,:)-x(j,:)
            d(1,:)=d(1,:)-L*nint(d(1,:)/L)
            r=d(1,1)**2+d(1,2)**2+d(1,3)**2
        if (r <= rc2) then
            do k=1,3
                ff(i,j,k)=48.0*e*sigma**12*d(1,k)*r**(-7)-24.0*e*&
		    sigma**6*d(1,k)*r**(-4)-48.0*e*sigma**12*d(1,k)*rc2**(-7)+24.0*e*sigma**6*d(1,k)*rc2**(-4)
                ff(j,i,k)=-ff(i,j,k)
            enddo
           E_V=E_V+4.0*e*(sigma**12/r**6-sigma**6/r**3) 
           
        else 
            do k=1,3
                ff(j,i,k)=0.0
                ff(i,j,k)=-ff(j,i,k)
            enddo
        endif
        enddo
    enddo
    

       
   
    Do k=1,3 
        Do i=1,npart
            Do j=1,npart    
                f_sum=f_sum+ff(i,j,k)
            enddo
            f(i,k)=f_sum
            f_sum=0.0
        enddo

    enddo
    
    E_V=E_V/npart
    
    
   
    
    return

    end subroutine

!> calculate kinetic energy 
    subroutine cal_ET(v,E_T)
    integer::i,j
    
    real(8),intent(in)::v(npart,3)
    real(8),intent(out)::E_T
    
    E_T=0.0
    
    Do i=1,npart
        Do j=1,3
            E_T=E_T+0.5*m*v(i,j)**2
        enddo
    enddo
    
    E_T=E_T/npart
       
    return

    end subroutine 
    
!> calculate velocity of particles 
    subroutine cal_V(f0,f,a,v0,v,E_T)
    integer::i,j
    real(8)::vv,vv1,vv2,p
    real(8),intent(in)::f0(npart,3),f(npart,3),v0(npart,3),a(1,3)
    real(8),intent(out)::v(npart,3),E_T
    
    E_T=0.0
    
    Do j=1,3
        Do i=1,npart
            v(i,j)=v0(i,j)+0.5*(f(i,j)+f0(i,j))*dt/m
            E_T=E_T+0.5*m*(v(i,j)**2)
        enddo
    enddo
    
    E_T=E_T/npart   
    
    p=freq*exp(-freq*dt)
    call random_seed()
    do i=1,npart
        call RANDOM_NUMBER(vv)
        if (vv <= p) then
        do j=1,3
            call RANDOM_NUMBER(vv1)
            call RANDOM_NUMBER(vv2)
            vv=sqrt(-2*log(vv1))*cos(2*pi*vv2)
            v(i,j)=vv*a(1,j)
        enddo
        endif
    enddo
    
    return

    end subroutine 

    end module

    
    
!***************************************************************
!> program main
program main
use constants
implicit none

real(8)::t,E_V,E_T
real(8)::x0(npart,3),v0(npart,3),f0(npart,3),x(npart,3),v(npart,3),f(npart,3),a(1,3)
integer::i
call inital(x0,v0,a)
call cal_ET(v0,E_T)
call cal_F(x0,f0,E_V)
open(unit=8,file='result.txt')
open(unit=9,file='velocity.txt')

t=0.0


do while (t.lt.tmax)
    x=x0+v0*dt+0.5/m*f0*dt**2
    
    call cal_F(x,f,E_V)
    call cal_V(f0,f,a,v0,v,E_T)
    
    if (MOD(t,(dt*nsample)) == 0) then  
    write (8,'(4(e16.8))')   t,E_T,E_V,E_T+E_V
    end if
    x0=x
    v0=v
    f0=f
    t=t+dt
    write(*,*) t
    
    if (t>=99E-13) then
    do i=1,npart
    write (9,'(3(e16.8))')   v(i,1),v(i,2),v(i,3)
    enddo
    endif


    enddo

close(8)
close(9)


end program main




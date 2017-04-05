!Write by Liao Yuanda
!
!###############################################################

module dou_ex
contains


    subroutine cau_Hf( L, N, M, Hf, t )
        !###############################################
        !(i,j) to (i, j_new); (i,j) to (i_new,j); N = i+(j-1)*L
        !###############################################
        implicit none
        integer :: L, N, M, i, j, i_new, j_new, N_ij, N_i1j, N_ij1
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) =0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                !period of boundary 
                i_new = mod(i+1, L)
                j_new = mod(j+1, L)

                N_ij = i*L+j+1
                N_i1j = i_new*L+j+1
                N_ij1 = i*L+j_new+1
                !Hf will be 0, if the bit is the same
                Hf(N_i1j, N_ij) = Hf(N_i1j, N_ij)-t
                Hf(N_ij1, N_ij) = Hf(N_ij1, N_ij)-t
                Hf(N_i1j+N, N_ij+N) = Hf(N_i1j+N, N_ij+N)-t
                Hf(N_ij1+N, N_ij+N) = Hf(N_ij1+N, N_ij+N)-t
            end do
        end do
    end subroutine cau_Hf

    subroutine cau_Hsf( L, N, M, Hsf, kesai, confi )
        implicit none
        integer :: L, N, M, i, j, k, N_ij, confi
        real*8 :: kesai, Hsf(M, M)
        do i = 1, M
            do j = 1, M
                Hsf(i,j) =0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                N_ij = i*L+j+1
                Hsf(N_ij, N_ij) = Hsf(N_ij, N_ij)-kesai*(ibits( confi, N_ij, 1 )-1/2) 
                Hsf(N_ij+N, N_ij+N) = Hsf(N_ij+N, N_ij+N)+kesai*(ibits( confi, N_ij, 1 )-1/2) 
            end do
        end do
        !do i = 1, M
        !    write(*,*) Hsf(i, :)
        !end do
    end subroutine cau_Hsf
    
    subroutine cau_H( M, Hf, Hsf, H )
        implicit none
        integer :: M, i, j
        real*8 Hf(M,M), Hsf(M,M), H(M, M)
        do i = 1, M
            do j = 1, M
                H(i,j) =0
            end do
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+Hsf(i,j)
            end do
        end do


    end subroutine cau_H
    subroutine cau_energy(M, H, energy )
        implicit none
        integer :: M, info, i, j
        real*8 :: H(M,M), energy
        real*8, dimension(M/2) :: Eig
        real*8, dimension(3*M/2) :: Work
        real*8, dimension(1, M/2) :: Vector
        real*8, dimension((M/2)*(M/2+1)/2) :: Hap
         do i = 1, M/2
             Eig(i) = 0
            do j = 1, M/2
                if (i <= j) then
                    Hap(int(i+(j-1)*j/2)) = H(i,j)
                end if
            end do
        end do
        call dspev('N', 'U', M/2, Hap, Eig, Vector, 1, Work, info)
        if (info /= 0) then
            write(*,*) info
        end if
        energy = 0.0
        do i = 1, M/2
            energy = energy + 2*Eig(i)
        end do
    end subroutine cau_energy
    
    subroutine cau_W(M, H, beta, Weight, Hap, Eig, Work, Vector )
        implicit none
        integer :: M, confi, info, i, j
        real*8 :: H(M,M), beta, Weight
        real*8 :: Hap(int(M*(M+1)/2)), Eig(M), Vector(1, M), Work(3*M)
         do i = 1, M
             Eig(i) = 0
            do j = 1, M
                if (i <= j) then
                    Hap(int(i+(j-1)*j/2)) = H(i,j)
                end if
            end do
        end do
        call dspev('N', 'U', M, Hap, Eig, Vector, 1, Work, info)
        if (info /= 0) then
            write(*,*) info
        end if
        Weight = 1.0
        !write(*,*) Eig(:)
        do i = 1, M
            Weight = Weight*( 1+exp( (-1)*beta*Eig(i) ) )
        end do
    end subroutine cau_W
    
    subroutine cau_Z(N, Weight, Z)

        implicit none
        integer :: N, i
        real*8 :: Weight(2**N), Z
        Z = 0.0
        do i = 1, 2**N
            Z = Z+Weight(i)
        end do
    end subroutine cau_Z

end module dou_ex

!subroutine init_random_seed()
!    integer :: i,n,clock
!    integer,dimension(:),allocatable :: seed
!    call random_seed(size=n)
!    allocate(seed(n))
!    call system_clock(count=clock)
!    seed=clock+37*(/(i-1,i=1,n)/)
!    call random_seed(PUT=seed) 
!    deallocate(seed)
!end subroutine init_random_seed
 

real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
    implicit none

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64

    ran64=ran64*mul64+add64
    ran=0.5d0+dmu64*dble(ran64)

end function ran
!----------------!

!---------------------!
subroutine initran(w)
!---------------------!
    implicit none

    integer(8) :: irmax
    integer(4) :: w,nb,b

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64
                        
    irmax=2_8**31
    irmax=2*(irmax**2-1)+1
    mul64=2862933555777941757_8
    add64=1013904243
    dmu64=0.5d0/dble(irmax)

    open(10,file='seed.in',status='old')

    read(10,*)ran64
    close(10)
    if (w.ne.0) then
        open(10,file='seed.in',status='unknown')
        write(10,*) abs((ran64*mul64)/5+5265361)
        close(10)
    end if

end subroutine initran
!----------------------!
program double_exchange
    use dou_ex 
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, kesai, the parameter if H
    !##########################################
    integer, parameter :: L = 4, N = 16, M = 32
    real*8, parameter :: t = 1
    integer :: i, j, k, ii, kk, num_MC, p, num_aver
    real*8 :: ran, s_single, s_aver, energy_aver, kesai, beta, kb, Tc, energy, x, Weight, Weight_new
    integer, allocatable, dimension(:) :: confi
    real*8, allocatable, dimension(:):: En
    real*8, dimension(M) :: Eig
    real*8, dimension(3*M) :: Work
    real*8, dimension(int(M*(M+1)/2)) :: Hap
    real*8, dimension(M,M) :: H, Hf, Hsf, Vector
    real*8, dimension(100) :: temp
    do i = 1, 100
        temp(i) = 0.05+0.01*i
    end do
    kesai = t
    
    call initran(1)
    open(unit=12, file='T-S.dat')
    num_MC = 6*2**N 
    allocate(confi(num_MC+1), En(num_MC+1))
do ii = 1, 100
    Tc = temp(ii)*t
    beta = 1/Tc 

    !write(*,*) 'b'
    call cau_Hf( L, N, M, Hf, t )
    !write(*,*) 'c'
    !do i = 1, M
    !    write(*,*) Hf(1,:)
    !end do
    x = ran()
    confi = 0
    En = 0
    confi(1) = int( x*(2**N))
    !write(*,*) 'd'
    i = 1
    do while (i <= num_MC)
    !write(*,*) 'e'
        call cau_Hsf( L, N, M, Hsf, kesai, confi(i))
        call cau_H( M, Hf, Hsf, H )
        call cau_energy(M, H, energy)
        call cau_W(M, H, beta, Weight, Hap, Eig, Work, Vector )
        !call random_seed
        x = ran()
        p = int(x*N) 
        !Write(*,*) p
        confi(i+1) = ibchng( confi(i), p )
        call cau_Hsf( L, N, M, Hsf, kesai, confi(i+1))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, Weight_new, Hap, Eig, Work, Vector )
        write(*,*) confi(i), confi(i+1)
        !write(*,*) exp((-1)*beta*(energy_new-energy))
        !call random_seed
        x = ran()
        write(*,*) x, '<?',  Weight_new/Weight, Weight
        Write(*,*) energy, (-1)/beta*log(Weight)
        if ( x <= Weight_new/Weight) then
            En(i) = energy
            i = i+1
            write(*,*) 'Step remainder:', num_MC-i
        end if

    !write(*,*) 'f'
    end do
    !open( unit = 11, file= 'confienergy.dat')
    !do j = 1, 200
    !    if (j < 150) then
    !        write(11,*)  j,  confi(int(j*num_MC/200)), En(int(j*num_mc/200))
    !    else
    !        Write(11,*)  j, confi(num_MC-200+j), En(num_MC-200+j)
    !    end if

    !end do
    !do i = 1, num_MC
    !    write(11,*) i, confi(i), En(i)
    !end do
    s_aver = 0
    energy_aver = 0
    num_aver = 2**N
    do j = 1, num_aver 
        s_single = 0
        do k = 1, N
            do kk = 1, N
                s_single = s_single+(ibits( confi(num_MC-j), k-1, 1 )-1/2)*2*(ibits( confi(num_MC-j), kk-1, 1 )-1/2)*2
            end do
        end do
        energy_aver = energy_aver + En(num_MC-j)
        S_single = s_single/N/N
        s_aver = s_aver+s_single
    end do
    s_aver = s_aver/num_aver
    energy_aver = energy_aver/num_aver
    Write(12,*) Tc, s_aver, energy_aver
    Write(*,*) 'step remainder:', 100-ii
    !pause 
end do

    deallocate(confi, En )
    write(*,*) 'end'



     

end program double_exchange

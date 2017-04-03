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

                N_ij = i*L+j
                N_i1j = i_new*L+j
                N_ij1 = i*L+j_new
                !Hf will be 0, if the bit is the same
                Hf(N_ij, N_i1j) = Hf(N_ij, N_i1j)-t
                Hf(N_ij, N_ij1) = Hf(N_ij, N_ij1)-t
                Hf(N_ij+N, N_i1j+N) = Hf(N_ij+N, N_i1j+N)-t
                Hf(N_ij+N, N_ij1+N) = Hf(N_ij+N, N_ij1+N)-t
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
                N_ij = i+j*L
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
    subroutine cau_W(M, H, beta, confi, Weight, Hap, Eig, Work, Vector )
        implicit none
        integer :: M, confi, info, i, j
        real*8 :: H(M,M), beta, Weight
        real*8 :: Hap(int(M*(M+1)/2)), Eig(M), Vector(1, M), Work(3*M)
         do i = 1, M
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



program double_exchange
    use dou_ex 
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, kesai, the parameter if H
    !##########################################
    integer :: L, N, i, j, k, num_MC, p, M
    real*8 :: t, kesai, beta, kb, Tc, energy, energy_new, x, Weight
    integer, allocatable, dimension(:) :: confi
    real*8, allocatable, dimension(:):: En
    real*8, allocatable, dimension(:) :: Hap, Eig, Work
    real*8, allocatable, dimension(:,:) :: H, Hf, Hsf, Vector
    L = 4 
    N = L*L
    M = 2*N
    t = 0.00001 
    kesai = 16.0*t
    Tc = 1*t
    beta = 1.0

    num_MC = 10*2**N 
    write(*,*) 'a'
    allocate( Hf(M, M), Hsf(M, M), H(M, M) )
    allocate(confi(num_MC+1), En(num_MC+1))
    allocate( Hap(int(M*(M+1)/2)), Eig(M), Vector(1,M), Work(3*M) )
    !write(*,*) 'b'
    call cau_Hf( L, N, M, Hf, t )
    !write(*,*) 'c'
    !do i = 1, M
    !    write(*,*) Hf(1,:)
    !end do
    call random_seed
    call random_number(x)
    confi = 0
    En = 0
    confi(1) = int( x*(2**N))
    !write(*,*) 'd'
    i = 1
    do while (i <= num_MC)
    !write(*,*) 'e'
        call cau_Hsf( L, N, M, Hsf, kesai, confi(i))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, confi(i), Weight, Hap, Eig, Work, Vector )
        write(*,*) Weight
        energy = (-1)*beta*log( Weight ) 
        call random_seed
        call random_number(x)
        p = int(x*N) 
        confi(i+1) = ibchng( confi(i), p )
        call cau_Hsf( L, N, M, Hsf, kesai, confi(i+1))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, confi(i+1), Weight, Hap, Eig, Work, Vector )
        energy_new = (-1)*beta*log( Weight ) 
        write(*,*) energy_new-energy, confi(i), confi(i+1)
        write(*,*) exp(-beta*(energy_new-energy))
        if ( 0.88 <= exp(-beta*(energy_new-energy))) then
            En(i) = energy
            i = i+1
            write(*,*) num_MC-i
        end if

    !write(*,*) 'f'
    end do
    open( unit = 11, file= 'confienergy.dat')
    do j = 1, 20
        write(11,*)  j,  confi(int(j*num_MC/20)), En(int(j*num_mc/20))
    end do

    write(*,*) 'Program end'
    deallocate(Eig)
    write(*,*) 'Eig end'
    deallocate(H)
    write(*,*) 'H end'
    deallocate(confi, En )
    write(*,*) 'end'


     
    write(*,*) 'Hf, Hsf end'

end program double_exchange

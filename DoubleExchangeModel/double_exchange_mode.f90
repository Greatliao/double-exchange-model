!Write by Liao Yuanda
!
!###############################################################

module dou_ex
contains

    subroutine init_basis( N, M, state)
        integer :: i, j, k, N, M, a, num_1, num_2, num_3, state(2**(3*N))
        M = 0
        do i = 0, N, 1
            do j = 0, N, 1
                do a = 0, 2**(3*N)-1, 1
                    num_1 = 0
                    num_2 = 0
                    num_3 = 0
                    do k = 0, N-1, 1
                        if ( btest(a, k+2*N) ) then
                            num_1 = num_1+1
                        end if
                        if ( btest(a, k+N) ) then
                            num_2 = num_2+1
                        end if
                        if ( btest(a, k) ) then
                            num_3 = num_3+1
                        end if
                    end do
                    if ( num_1 == i .AND. num_2 == N-i .AND. num_3 == j) then
                        M = M+1
                        state(M) = a
                    end if
                end do
            end do
        end do
        write(*,*) 2**(3*N)
        Write(*,*) 2**(2*N)
        write(*,*) M
    end subroutine init_basis

    subroutine cau_Hf( L, N, M, basis, Hf, t )
        !###############################################
        !(i,j) to (i, j_new); (i,j) to (i_new,j); N = i+(j-1)*L
        !###############################################
        implicit none
        integer :: L, N, M, basis(M), i, j, k, i_new, j_new, N_ij, N_i1j, N_ij1, b, e
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) =0
            end do
        end do

        do k = 1, M
            do i = 0, L-1, 1
                do j = 0, L-1, 1
                    !period of boundary 
                    i_new = mod(i+1, L)
                    j_new = mod(j+1, L)

                    N_ij = i+j*L
                    N_i1j = i_new+j*L
                    N_ij1 = i+j_new*L
                    !Hf will be 0, if the bit is the same
                    if ( ( ibits( basis(k), 2*N+N_ij, 1 ) /= ibits( basis(k), 2*N+N_i1j, 1 ) ) .AND. &
                        ( ibits( basis(k), N+N_ij, 1 ) /= ibits( basis(k), N+N_i1j, 1 ) ) ) then
                        b = ibchng( basis(k), 2*N+N_ij )
                        b = ibchng( b, 2*N+N_i1j)
                        b = ibchng( b, N+N_ij)
                        b = ibchng( b, N+N_i1j)
                        do e = 1, M
                            if ( b == basis(e) ) then
                                Hf(e,k) = Hf(e,k)-t
                            end if
                            exit
                        end do
                    end if
                    if ( ( ibits( basis(k), 2*N+N_ij, 1 ) /= ibits( basis(k), 2*N+N_ij1, 1 ) ) .AND. &
                        ( ibits( basis(k), N+N_ij, 1 ) /= ibits( basis(k), N+N_ij1, 1 ) ) ) then
                        b = ibchng( basis(k), 2*N+N_ij )
                        b = ibchng( b, 2*N+N_ij1)
                        b = ibchng( b, N+N_ij)
                        b = ibchng( b, N+N_ij1)
                        do e = 1, M
                            if ( b == basis(e) ) then
                                Hf(e,k) = Hf(e,k)-t
                            end if
                            exit
                        end do
                    end if
                end do
            end do
        end do
    end subroutine cau_Hf

    subroutine cau_Hs( L, M, basis, Hs, J_para)
        implicit none
        integer :: L, M, basis(M), i, j, k, i_new, j_new, N_ij, N_i1j, N_ij1
        real*8 :: J_para, Hs(M, M)
        do i = 1, M
            do j = 1, M
                Hs(i,j) =0
            end do
        end do
        do k = 1, M
            do i = 0, L-1, 1
                do j = 0, L-1, 1
                    !period of boundary 
                    i_new = mod(i+1, L)
                    j_new = mod(j+1, L)

                    N_ij = i+j*L
                    N_i1j = i_new+j*L
                    N_ij1 = i+j_new*L
                    
                    Hs(k,k) = Hs(k,k)-J_para*(ibits( basis(k), N_ij, 1 )-1/2)*(ibits( basis(k), N_i1j, 1 )-1/2) 
                    Hs(k,k) = Hs(k,k)-J_para*(ibits( basis(k), N_ij, 1 )-1/2)*(ibits( basis(k), N_ij1, 1 )-1/2) 
                end do
            end do
        end do
    end subroutine cau_Hs

    subroutine cau_Hsf( L, N, M, basis, Hsf, kesai )
        implicit none
        integer :: L, N, M, basis(M), i, j, k, N_ij
        real*8 :: kesai, Hsf(M, M)
        do i = 1, M
            do j = 1, M
                Hsf(i,j) =0
            end do
        end do

        do k = 1, M
            do i = 0, L-1, 1
                do j = 0, L-1, 1
                    N_ij = i+j*L

                    Hsf(k,k) = Hsf(k,k)-kesai*(ibits( basis(k), N_ij, 1 )-1/2)*(ibits( basis(k), 2*N+N_ij, 1 )-ibits( basis(k), N+N_ij, 1 ) ) 
                end do
            end do
        end do
    end subroutine cau_Hsf
    
    subroutine cau_H( M, Hf, Hs, Hsf, H )
        implicit none
        integer :: M, i, j
        real*8 Hf(M,M), Hs(M,M), Hsf(M,M), H(M, M)
        do i = 1, M
            do j = 1, M
                H(i,j) =0
            end do
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+Hs(i,j)+Hsf(i,j)
            end do
        end do

    end subroutine cau_H

    

end module dou_ex

program double_exchange

    use dou_ex
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !M, size of basis
    !t, J_para, kesai, the parameter if H
    !##########################################
    integer :: L, N, M, i, j, info
    real*8 :: t, J_para, kesai
    integer, allocatable, dimension(:) :: basis, state
    real*8, allocatable, dimension(:) :: Hap, Eig, Work
    real*8, allocatable, dimension(:,:) :: H, Hf, Hs, Hsf, Vector
    open(unit = 12, file = 'H.dat', status = 'replace', action = 'write')
    open(unit = 11, file = 'basis.dat',  status = 'replace', action = 'write')
    open(unit = 13, file = 'Eig.bat',  status = 'replace', action = 'write')
    L = 2
    N = L*L
    M = 0
    t = 1.0
    J_para = 1.0
    kesai = 1.0
    !init the basis
    allocate(state(2**(3*N)))
    call init_basis( N, M, state)
    allocate( basis(M), Hs(M, M), Hf(M, M), Hsf(M, M), H(M, M) )

    write(*,*) 'init_basis'
    do i = 1, M
        basis(i) = state(i)
        Write(11,*) basis(i)
    end do
    deallocate(state)

    !caculate the H
    call cau_Hf( L, N, M, basis, Hf, t )
    write(*,*) 'cau_Hf'
    call cau_Hs( L, M, basis, Hs, J_para ) 
    write(*,*) 'cau_Hs'
    call cau_Hsf( L, N, M, basis, Hsf, kesai ) 
    write(*,*) 'cau_Hsf'
    call cau_H( M, Hf, Hs, Hsf, H ) 
    deallocate(Hf, Hs, Hsf)
    
    write(*,*) 'cau_H'
    do i = 1, M
        write(12,*) H(i,:)
    end do
    allocate( Hap(int(M*(M+1)/2)), Eig(M), Vector(1,M), Work(3*M) )
    Hap = 0
    Eig = 0
    do i = 1, M
        do j = 1, M
            Vector(i,j) = 0
        end do
    end do

    do i = 1, M
        do j = 1, M
            if (i <= j) then
                Hap(int(i+(j-1)*j/2)) = H(i,j)
            end if
        end do
    end do

    call dspev('N', 'U', M, Hap, Eig, Vector, M, Work, info)
    write(*,*) 'info = ', info
    do i = 1, M
        write(13,*) Eig(i)
    end do
    close(12)
    close(13)
    close(11)
    deallocate(basis)
    deallocate(Eig)
    deallocate(H)
end program double_exchange

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
        integer :: L, N, M, i, j, k, i_minus, i_plus, j_minus, j_plus, k_minus, k_plus, N_ij, N_i_plus, N_i_minus, N_j_minus, N_j_plus, N_k_minus, N_k_plus
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                do k = 0, L-1, 1
                !period of boundary 
                i_plus = mod(i+1, L)
                j_plus = mod(j+1, L)
                k_plus = mod(k+1, L)
                if (i==0) then
                    i_minus = L-1
                else
                    i_minus = i-1
                end if 
                if (j==0) then
                    j_minus = L-1
                else
                    j_minus = j-1
                end if 
                if (k==0) then
                    k_minus = L-1
                else
                    k_minus = j-1
                end if
                N_ij = i*L*L+j*L+k+1
                N_i_minus = i_minus*L*L+j*L+k+1
                N_i_plus = i_plus*L*L+j*L+k+1
                N_j_minus = i*L*L+j_minus*L+k+1
                N_j_plus = i*L*L+j_plus*L*k+1
                N_k_minus = i*L*L+j*L+k_minus+1
                N_k_plus = i*L*L+j*L+k_plus+1

                !Hf will be 0, if the bit is the same
                Hf(N_ij, N_i_minus) = Hf(N_ij, N_i_minus)-t
                Hf(N_ij, N_i_plus) = Hf(N_ij, N_i_plus)-t
                Hf(N_ij, N_j_minus) = Hf(N_ij, N_j_minus)-t
                Hf(N_ij, N_j_plus) = Hf(N_ij, N_j_plus)-t
                Hf(N_ij, N_k_minus) = Hf(N_ij, N_k_minus)-t
                Hf(N_ij, N_k_plus) = Hf(N_ij, N_k_plus)-t

                Hf(N_ij+N, N_i_minus+N) = Hf(N_ij+N, N_i_minus+N)-t
                Hf(N_ij+N, N_i_plus+N) = Hf(N_ij+N, N_i_plus+N)-t
                Hf(N_ij+N, N_j_minus+N) = Hf(N_ij+N, N_j_minus+N)-t
                Hf(N_ij+N, N_j_plus+N) = Hf(N_ij+N, N_j_plus+N)-t
                Hf(N_ij+N, N_k_minus+N) = Hf(N_ij+N, N_k_minus+N)-t
                Hf(N_ij+N, N_k_plus+N) = Hf(N_ij+N, N_k_plus+N)-t
                end do
            end do
        end do
    end subroutine cau_Hf

    subroutine cau_Hsf( N, M, Hsf, J_para, fai, theta )
        implicit none
        integer :: N, M, i, j
        real*8 :: J_para, fai(N), theta(N)
        complex*16 :: Hsf(M, M)
        do i = 1, M
            do j = 1, M
                Hsf(i,j) = (0.0, 0.0)
            end do
        end do

        do i = 1, N 
            Hsf(i, i) = Hsf(i, i)-J_para*cos(theta(i)) 
            Hsf(i+N, i+N) = Hsf(i+N, i+N)+J_para*cos(theta(i)) 
            Hsf(i+N, i) = Hsf(i+N, i)-J_para*sin(theta(i))*exp( cmplx(0.0, fai(i)) ) 
            Hsf(i, i+N) = Hsf(i, i+N)-J_para*sin(theta(i))*exp( cmplx(0.0, fai(i)*(-1)) ) 

            !write(*,*) (ibits( confi, i-1, 1 )*1.0-1.0/2) 
            !write(*,*) Hsf(i, i)
        end do
        !write(*,*) fai(:)
        !write(*,*) theta(:)
        !write(*,*) "end"
        !do i = 1, M
        !    write(*,*) Hsf(i, :)
        !end do
    end subroutine cau_Hsf
    
    subroutine cau_H( M, Hf, Hsf, H )
        implicit none
        integer :: M, i, j
        real*8 :: Hf(M,M)
        complex*16 :: Hsf(M,M), H(M, M)
        do i = 1, M
            do j = 1, M
                H(i,j) =0.0
            end do
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+Hsf(i,j)
            end do
        end do


    end subroutine cau_H
    
    subroutine cau_W(M, H, beta, ln_Weight, Eig, Work, Rwork, miu )
        implicit none
        integer :: M, info, i, j
        real*8 :: beta, ln_Weight, miu
        complex*16 :: H(M,M), Work(3*M)
        real*8 :: Eig(M),  Rwork(3*M-2)

        call zheev('N', 'U', M, H, M, Eig, Work, 3*M, Rwork, info)
        if (info /= 0) then
            write(*,*) 'info_erro =', info
        end if
        ln_Weight = 0.0 
        !write(*,*) Eig(:)
        do i = 1, M
            ln_Weight = ln_Weight+log( 1+exp( (-1)*beta*(Eig(i)-miu) ) )
        end do
    end subroutine cau_W

end module dou_ex

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

!function to less dimension
integer function di(d)
    integer :: d
    if (d == 1) then
        di = 2
    end if
    if ( d == 2) then
        di =1
    end if
end function

program double_exchange
    use dou_ex 
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_para, the parameter if H
    !##########################################
    integer, parameter :: L = 4, N = 64, M = 128, Num_T = 100, num_MC = 10**6
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, p, num_aver, d, di, dd
    real*8 :: sx, sy, sz, ran, J_para, beta, miu, Tc, x, ln_Weight, ln_Weight_new_fai, ln_Weight_new_theta
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(num_MC+1) :: En, ss, s_M
    real*8, dimension(M) :: Eig
    real*8, dimension(3*M-2) :: Rwork
    real*8, dimension(M,M) :: Hf
    complex*16, dimension(3*M) :: Work
    complex*16, dimension(M,M) :: H, Hsf 
    J_para = 8*t
    miu = -8*t
    d = 1

    
    call initran(1)
    call cau_Hf( L, N, M, Hf, t )

    open(unit=12, file='s_all.dat')
    open(unit=13, file='En_all.dat')
    open(unit=14, file='s.dat')
    open(unit=15, file='En.dat')
    Tc = 0.5*t
    beta = 1/Tc

    !write(*,*) 'b'
    !write(*,*) 'c'
    !do i = 1, M
    !    write(*,*) Hf(1,:)
    !end do
    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI 
    end do
    !write(*,*) 'd'
    do i = 1, num_MC 
    !write(*,*) 'e'
        dd = di(d)
        call cau_Hsf( N, M, Hsf, J_para, fai(d,:), theta(d,:))
        call cau_H( M, Hf, Hsf, H)
        call cau_W(M, H, beta, ln_Weight, Eig, Work, Rwork , miu )
        En(i) = (-1)*ln_Weight/beta 
        !call random_seed
        x = ran()
        p = int(x*N)+1 
        !Write(*,*) p
        fai(dd,:) = fai(d,:)
        theta(dd,:) = theta(d,:)
        x = ran()
        fai(dd, p) = fai(d,p)+(x-0.5)*2*2*PI/12
        if ( fai(dd,p) < 0 ) then
            fai(dd,p) = fai(dd,p)+2*PI
        end if
        if ( fai(dd,p) > 2*PI ) then
            fai(dd,p) = fai(dd,p)-2*PI
        end if
        !call cau_Hsf( N, M, Hsf, J_para, fai(dd,:), theta(dd,:))
        !call cau_H( M, Hf, Hsf, H, miu )
        !call cau_W(M, H, beta, ln_Weight_new_fai, Eig, Work, Rwork )

        !x = ran()
        !if ( x > exp( ln_Weight_new_fai-ln_Weight ) ) then
        !    fai(dd, p) = fai(d, p)
        !    ln_Weight_new_fai = ln_Weight
        !end if

        x = ran()
        theta(dd, p) = theta(d,p)+(x-0.5)*2*PI/6
        if ( theta(dd,p) < 0 ) then
            theta(dd,p) = theta(dd,p)+PI
        end if
        if (theta(dd,p) > PI ) then
            theta(dd,p) = theta(dd,p)-PI
        end if
        call cau_Hsf( N, M, Hsf, J_para, fai(dd,:), theta(dd,:))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, ln_Weight_new_theta, Eig, Work, Rwork , miu)

        x = ran()
        !write(*,*) ln_Weight, ln_Weight_new_theta
        !write(*,*) x, '<?',  exp(ln_Weight_new_theta-ln_Weight)
        !Write(*,*) energy, (-1)/beta*log(ln_Weight)
        if ( x > exp( ln_Weight_new_theta-ln_Weight ) ) then
            fai(dd,p) = fai(d,p)
            theta(dd, p) = theta(d, p)
        end if
        ss(i) = 0
        sx = 0
        sy = 0 
        sz = 0
        s_M(i) = 0
        do j = 1, N
            sx = sx+sin(theta(d,j))*cos(fai(d,j))
            sy = sy+sin(theta(d,j))*sin(fai(d,j))
            sz = sz+cos(theta(d,j))
            do k = 1, N
                ss(i) = ss(i)+sin(theta(d,j))*cos(fai(d,j))*sin(theta(d,k))*cos(fai(d,k))
                ss(i) = ss(i)+sin(theta(d,j))*sin(fai(d,j))*sin(theta(d,k))*sin(fai(d,k))
                ss(i) = ss(i)+cos(theta(d,j))*cos(theta(d,k))
            end do
        end do
        s_M(i) = sqrt(sx**2+sy**2+sz**2)/N/N
        ss(i) = ss(i)/N
        write(12,*) i, ss(i), s_M(i)
        write(13,*) i, En(i)
        d = dd
        write(*,*) i/(num_MC/100)
    end do
    do i = 1, 10000
        write(14,*) i, ss(num_MC-10000+i), s_M(num_MC-10000+i)
        write(15,*) i, En(num_MC-10000+i)
    end do
    write(*,*) 'end'
end program double_exchange

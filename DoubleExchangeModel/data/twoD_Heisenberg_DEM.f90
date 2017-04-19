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
        integer :: L, N, M, i, j, i_left, i_right, j_down, j_up, N_ij, N_i_right, N_i_left, N_j_down, N_j_up 
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                !period of boundary 
                i_right = mod(i+1, L)
                j_down = mod(j+1, L)
                if (i==0) then
                    i_left = L-1
                else
                    i_left = i-1
                end if
                if (j==0) then
                    j_up = L-1
                else
                    j_up = j-1
                end if

                N_ij = i*L+j+1
                N_i_right = i_right*L+j+1
                N_i_left = i_left*L+j+1
                N_j_down = i*L+j_down+1
                N_j_up = i*L+j_up+1
                !Hf will be 0, if the bit is the same
                Hf(N_ij, N_i_right) = Hf(N_ij, N_i_right)-t
                Hf(N_ij, N_i_left) = Hf(N_ij, N_i_left)-t
                Hf(N_ij, N_j_down) = Hf(N_ij, N_j_down)-t
                Hf(N_ij, N_j_up) = Hf(N_ij, N_j_up)-t

                Hf(N_ij+N, N_i_right+N) = Hf(N_ij+N, N_i_right+N)-t
                Hf(N_ij+N, N_i_left+N) = Hf(N_ij+N, N_i_left+N)-t
                Hf(N_ij+N, N_j_down+N) = Hf(N_ij+N, N_j_down+N)-t
                Hf(N_ij+N, N_j_up+N) = Hf(N_ij+N, N_j_up+N)-t
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
                H(i,j) =(0.0, 0.0)
            end do
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = H(i,j)+Hf(i,j)+Hsf(i,j)
            end do
        end do
        !do i = 1, M
        !    write(*,*) 'comn**********', i 
        !write(*,*) H(i,:)
        !end do
    end subroutine cau_H
    
    subroutine cau_W(M, H, beta, ln_Weight, Eig, Work, Rwork, miu )
        implicit none
        integer :: M, info, i, j
        real*8 :: beta, ln_Weight, miu
        complex*16 :: H(M,M), Work(3*M)
        real*8 :: Eig(M),  Rwork(3*M-2)
        
        call zheev('N', 'U', M, H, M, Eig, Work, 3*M, Rwork, info)
        if (info /= 0) then
            write(*,*) 'info_erro = ', info
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
    integer, parameter :: L = 4, N = 16, M = 32, num_MC = 10**6, num_T = 50
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, p, d, di, dd, ii, num_aver
    real*8 :: ran, J_para, beta, miu, kb, Tc, x, sx, sy, sz, ln_Weight,ln_Weight_test,  ln_Weight_new_fai, ln_Weight_new_theta
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(num_MC) :: En, ss, s_M
    real*8, dimension(M) :: Eig
    real*8, dimension(3*M-2) :: Rwork
    real*8, dimension(M,M) :: Hf
    complex*16, dimension(3*M) :: Work
    complex*16, dimension(M,M) :: H, Hsf 
    real*8, dimension(num_T) :: temp
    J_para = 8*t
    miu = -8*t
    d = 1
    do i = 1, num_T
        temp(i) = 0.01+0.01*i
    end do
    call initran(1)
    open(unit=14,file='T-S.dat')
    call cau_Hf( L, N, M, Hf, t )

do ii = 1, num_T
    Tc = temp(ii)*t
    beta = 1/Tc

    do i = 1, N
        x = ran()
        fai(1,i) = x*2*PI
        x = ran()
        theta(1,i) = x*PI 
    end do
    do i = 1, num_MC 
        dd = di(d)
        call cau_Hsf( N, M, Hsf, J_para, fai(d,:), theta(d,:))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, ln_Weight, Eig, Work, Rwork, miu )
        !write(*,*) '######################################'
        !call cau_Hsf( N, M, Hsf, J_para, fai(d,:), theta(d,:))
        !call cau_H( M, Hf, Hsf, H, miu )
        !call cau_W(M, H, beta, ln_Weight_test, Eig, Work, Rwork )
        !write(*,*) ln_Weight, ln_Weight_test 
        !pause
        En(i) = (-1)*ln_Weight/beta
        x = ran()
        p = int(x*N)+1 
        !write(*,*) p
        fai(dd,:) = fai(d,:)
        theta(dd,:) = theta(d,:)
        x = ran()
        fai(dd, p) = fai(dd,p)+(x-0.5)*2*2*PI/12
        if ( fai(dd,p) < 0 ) then
            fai(dd,p) = fai(dd,p)+2*PI
        end if
        if ( fai(dd,p) > 2*PI ) then
            fai(dd,p) = fai(dd,p)-2*PI
        end if
     
        !call cau_Hsf( N, M, Hsf, J_para, fai(dd,:), theta(dd,:))
        !call cau_H( M, Hf, Hsf, H )
        !call cau_W(M, H, beta, ln_Weight_new_fai, Eig, Work, Rwork, miu )

        !x = ran()
        !write(*,*) 'fai=', fai(i,:)
        !write(*,*) x, "<?", exp(ln_Weight_new_fai-ln_Weight)
        !if ( x > exp( ln_Weight_new_fai-ln_Weight ) ) then
        !    fai(dd, p) = fai(d, p)
        !    ln_Weight_new_fai = ln_Weight
        !end if
        !write(*,*) ln_Weight, ln_Weight_new_fai

        x = ran()
        theta(dd, p) = theta(dd,p)+(x-0.5)*2*PI/6
        if ( theta(dd,p) < 0 ) then
            theta(dd,p) = theta(dd,p)+PI
        end if
        if ( theta(dd,p) > PI ) then
            theta(dd,p) = theta(dd,p)-PI
        end if

        call cau_Hsf( N, M, Hsf, J_para, fai(dd,:), theta(dd,:))
        call cau_H( M, Hf, Hsf, H )
        call cau_W(M, H, beta, ln_Weight_new_theta, Eig, Work, Rwork, miu )

        x = ran()
        !write(*,*) ln_Weight, ln_Weight_new_theta
        !write(*,*) 'theta=', theta(i,:)
        !write(*,*) x, "<?", exp(ln_Weight_new_theta-ln_Weight_new_fai)
        !write(*,*) p
        !write(*,*) fai(d,:)
        !write(*,*) fai(dd,:)
        !write(*,*) theta(d,:)
        !write(*,*) theta(dd,:)
        if ( x > exp( ln_Weight_new_theta-ln_Weight ) ) then
            fai(dd,p) = fai(d,p)
            theta(dd, p) = theta(d, p)
        end if
        ss(i) = 0
        s_M(i) = 0
        sx = 0
        sy = 0
        sz = 0
        do j = 1, N
            sx = sx+sin(theta(d,j))*cos(fai(d,j))
            sy = sy+sin(theta(d,j))*sin(fai(d,j))
            sz = sz+cos(theta(d,j))
            do k = 1, N
                ss(i)=ss(i)+sin(theta(d,j))*cos(fai(d,j))*sin(theta(d,k))*cos(fai(d,k))
                ss(i)=ss(i)+sin(theta(d,j))*sin(fai(d,j))*sin(theta(d,k))*sin(fai(d,k))
                ss(i)=ss(i)+cos(theta(d,j))*cos(theta(d,k))
            end do
        end do
        s_M(i) = sqrt(sx**2+sy**2+sz**2)/N
        ss(i) = ss(i)/N/N
        !En(i) = energy 
        !write(*,*) i/(num_MC/100)
        !write(*,*) d
        d = dd
    end do
    num_aver = 10000
    sx = 0
    sy = 0
    do i = 1, num_aver
       sx = sx+ss(num_MC-num_aver+i)
       sy = sy+s_M(num_MC-num_aver+i)
    end do
    sx = sx/num_aver
    sy = sy/num_aver
    write(14,*) ii, sx, sy
    write(*,*) num_T-ii
end do
    
    write(*,*) 'end'
end program double_exchange

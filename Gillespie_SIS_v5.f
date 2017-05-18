      !====================================================================
      !                        Walter Ullon - 08/2016
      !                           Gillespie SIS v5
      ! This program is the fortran version of Gillespie_SIS_v3 ported from
      ! MATLAB, without the implementation of control measures. 
      ! A Montecarlo simulation of the SIS epidemic model
      ! This version saves the extinction time for a given number of runs,
      ! and returns the mean extinction time (MTE), given R0.
      !====================================================================

      program Gillespie_SIS

      implicit none

      integer*8, parameter :: i_max = 10.0**12   ! max. number of iterations
      integer*8, parameter :: u_max = 125        ! max. number of iterations
      integer*8 :: t_ext, I_init, S_init, S_now, I_now, i, u, j
      integer :: N, endemic_I, endemic_S
      real*8 :: t, a0, tau, st, r1, r2

      real, parameter :: beta = 1000           ! contact rate
      real, parameter :: gammma = 99.98        ! recovery rate
      real, parameter :: mu = 0.02             ! birth/death rate
      real, parameter :: R0 = beta/(gammma+mu) ! reproductive number
     
      real, dimension(u_max) :: t_new  ! time array
      real, dimension(5) :: events     ! events array

      !========================================================
      ! RNG: initialize, seed using system-time - DO NOT TOUCH!
      !========================================================
      integer :: values(1:8), k
      integer, dimension(:), allocatable :: seed
      real(8) :: r
      call date_and_time(values=values)
      call random_seed(size=k)
      allocate(seed(1:k))
      seed(:) = values(8)
      call random_seed(put=seed)

      open(unit=3, file='run8.dat') ! data storage file

      N = 40                    ! initial population size
      j = 1                     ! time array index

      endemic_S = nint(N/R0)          ! endemic steady state (S)
      endemic_I = nint(N*(1-(1/R0)))  ! endemic steady state (I)
      
      !=============================================
      ! main loop
      !=============================================
      do u=1, u_max
         t = 0
         t_ext = 0              ! extiction time
         
         I_init = 2             ! initial infected pop.
         S_init = N-I_init      ! initial susceptible pop.

         S_now = S_init            ! holding variable for susceptible
         I_now = I_init            ! holding variable for infected

         !===========================================
         ! sub loop
         !===========================================
         do i=1, i_max
            call random_number(r1)       ! assign randon number to r1
            call random_number(r2)       ! assign random number to r2

         !==========================================
         ! generate events
         !==========================================
            events(1) = mu*N                 ! Birth(S)
            events(2) = mu*S_now             ! Death(S)
            events(3) = mu*I_now             ! Death(I)
            events(4) = (beta*S_now*I_now)/N ! Infection
            events(5) = gammma*I_now         ! Recovery

            a0 = events(1)+events(2)+events(3)+events(4)+events(5)

            st = r2*a0             ! determine which event occurs
            tau = LOG(1/r1)*(1/a0) ! determine when said event occurs
            t = t + tau            ! update time

         !============================================
         ! update the populations
         !============================================
            if (st .le. events(1)) then
                S_now = S_now + 1
            else if (st .gt. events(1) .AND. st .le. 
     #(events(1) + events(2)))  then
                S_now = S_now - 1
            else if (st .gt. (events(1) + events(2)) .AND. st .le. 
     #(events(1) + events(2) + events(3))) then
                I_now = I_now - 1
            else if (st .gt. (events(1) + events(2) + events(3)) .AND. 
     #st .le. (events(1) + events(2) + events(3) + events(4))) then
               S_now = S_now - 1
               I_now = I_now + 1
            else
               S_now = S_now + 1
               I_now = I_now - 1
            end if
         
         !============================================
         ! capture extinction & write to file
         !============================================
            if(I_now .le. 0) then  ! if extinct, terminate
              ! print *, j
               t_new(j) = t
               write(3,*) t_new(j)
               j = j+1
               exit
            end if 

         end do                 ! subloop
      
      end do                    ! main loop

      print*, " "
      print*, "Pop. size = ", N
      print*, "R0 = ", R0
      print*, "endemic I = ", endemic_I
      print*, "endemic S = ", endemic_S
      !print*, "MTE = ", (log(sum(t_new)/u_max)/N)
      print*, "MTE = ", (sum(t_new)/j)
      print*, "Number of extinctions captured = ", size(t_new)
      print*, " "

       end program Gillespie_SIS

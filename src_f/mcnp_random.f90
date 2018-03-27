module mcnp_random
implicit none
  !=======================================================================
  ! Description:
  !  mcnp_random.F90 -- random number generation routines
  !=======================================================================
  ! Simplified, only 1 63 bit generator (#2)
  !=======================================================================
  !  This module contains:
  !
  !   * Constants for the RN generator, including initial RN seed for the
  !     problem & the current RN seed
  !
  !   * MCNP interface routines:
  !     - random number function:           rang()
  !     - RN initialization for problem:    RN_init_problem
  !     - RN initialization for particle:   RN_init_particle
  !     - get info on RN parameters:        RN_query
  !     - get RN seed for n-th history:     RN_query_first
  !     - set new RN parameters:            RN_set
  !     - skip-ahead in the RN sequence:    RN_skip_ahead
  !     - Unit tests:        RN_test_basic, RN_test_skip, RN_test_mixed
  !
  !   * For interfacing with the rest of MCNP, arguments to/from these
  !     routines will have types of I8 or I4.
  !     Any args which are to hold random seeds, multipliers,
  !     skip-distance will be type I8, so that 63 bits can be held without
  !     truncation.
  !
  ! Revisions:
  ! * 10-04-2001 - F Brown, initial mcnp version
  ! * 06-06-2002 - F Brown, mods for extended generators
  ! * 12-21-2004 - F Brown, added 3 of LeCuyer's 63-bit mult. RNGs
  ! * 01-29-2005 - J Sweezy, Modify to use mcnp modules prior to automatic
  !                io unit numbers.
  ! * 12-02-2005 - F Brown, mods for consistency with C version
  ! * 12-12-2017 - F Brown, update to generator #2
  !=======================================================================

  !-------------------
  ! MCNP output units
  !-------------------
  integer, parameter ::   iuo = 6
  integer, parameter ::  jtty = 6

  PRIVATE
  !---------------------------------------------------
  ! Kinds for LONG INTEGERS (64-bit) & REAL*8 (64-bit)
  !---------------------------------------------------
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: I8 = selected_int_kind(18)

  !-----------------------------------
  ! Public functions and subroutines for this module
  !-----------------------------------
  PUBLIC ::  rang
  PUBLIC ::  RN_init_problem
  PUBLIC ::  RN_init_particle
  PUBLIC ::  RN_test_basic

  !-----------------------------------------------------------------
  !   * Linear multiplicative congruential RN algorithm:
  !
  !            RN_SEED = RN_SEED*RN_MULT + RN_ADD  mod RN_MOD
  !
  !   * Default values listed below will be used, unless overridden
  !-----------------------------------------------------------------
  integer(I8), SAVE      :: RN_SEED0  = 1_I8
  integer(I8), parameter :: RN_MULT   = 9219741426499971445_I8
  integer(I8), parameter :: RN_ADD    = 1_I8
  integer,     parameter :: RN_BITS   = 63
  integer(I8), parameter :: RN_STRIDE = 152917_I8
  integer(I8), parameter :: RN_PERIOD = ishft( 1_I8, RN_BITS )
  integer(I8), parameter :: RN_MOD    = ishft( 1_I8, RN_BITS )
  integer(I8), parameter :: RN_MASK   = ishft( not(0_I8), RN_BITS-64 )
  integer(I8), parameter :: RN_SHIFT  = 53 - RN_BITS
  real(R8),    parameter :: RN_NORM   = 2._R8**(-53)

  !------------------------------------
  ! Private data for a single particle
  !------------------------------------
  integer(I8) :: RN_SEED   = 1_I8 ! current seed

  common                /RN_THREAD/   RN_SEED
  save                  /RN_THREAD/
  !$OMP THREADprivate ( /RN_THREAD/ )

  !---------------------------------------------------------------------
  ! Reference data:  Seeds for case of init.seed = 1,
  !                  Seed numbers for index 1-5, 123456-123460
  !---------------------------------------------------------------------
  integer(I8), parameter, dimension(10) ::  RN_CHECK = [ &
    ! ***** 2 *****
    &  9219741426499971446_I8,  666764808255707375_I8, 4935109208453540924_I8, &
    &  7076815037777023853_I8, 5594070487082964434_I8, 7069484152921594561_I8, &
    &  8424485724631982902_I8,   19322398608391599_I8, 8639759691969673212_I8, &
    &  8181315819375227437_I8 ]
  !---------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------

  function rang()
    ! MCNP random number generator
    implicit none
    real(R8) ::  rang

    RN_SEED  = iand( iand(RN_MULT*RN_SEED,RN_MASK) + RN_ADD, RN_MASK) 
    rang     = max( ishft(RN_SEED,RN_SHIFT), 1_I8) * RN_NORM

    return
  end function rang

  !-------------------------------------------------------------------

  function RN_skip_ahead( seed, skip )
    ! advance the seed "skip" RNs:   seed*RN_MULT^n mod RN_MOD
    implicit none
    integer(I8) :: RN_skip_ahead
    integer(I8), intent(in)  :: seed, skip
    integer(I8) :: nskip, gen, g, inc, c, gp, rn, seed_old

    seed_old = seed
    ! add period till nskip>0
    nskip = skip
    do while( nskip<0_I8 )
      if( RN_PERIOD>0_I8 ) then
        nskip = nskip + RN_PERIOD
      else
        nskip = nskip + RN_MASK
        nskip = nskip + 1_I8
      endif
    enddo

    ! get gen=RN_MULT^n,  in log2(n) ops, not n ops !
    nskip = iand( nskip, RN_MASK )
    gen   = 1
    g     = RN_MULT
    inc   = 0
    c     = RN_ADD
    do while( nskip>0_I8 )
      if( btest(nskip,0) )  then
        gen = iand( gen*g, RN_MASK )
        inc = iand( inc*g, RN_MASK )
        inc = iand( inc+c, RN_MASK )
      endif
      gp    = iand( g+1,  RN_MASK )
      g     = iand( g*g,  RN_MASK )
      c     = iand( gp*c, RN_MASK )
      nskip = ishft( nskip, -1 )
    enddo
    rn = iand( gen*seed_old, RN_MASK )
    rn = iand( rn + inc, RN_MASK )
    RN_skip_ahead = rn
    return
  end function RN_skip_ahead

  !-------------------------------------------------------------------

  subroutine RN_init_problem( new_seed,  print_info )
    ! * initialize MCNP random number parameters for problem,
    !   based on user input.  This routine should be called
    !   only from the main thread, if OMP threading is being used.
    !
    ! * check on size of long-ints & long-int arithmetic
    ! * check the multiplier
    ! * advance the base seed for the problem
    ! * set the initial particle seed
    implicit none
    integer(I8),           intent(in) :: new_seed
    integer,     optional, intent(in) :: print_info
    character(len=20) :: printseed
    integer(I8)       ::  itemp1, itemp2, itemp3, itemp4

    ! set defaults, override if input supplied: seed
    if( new_seed>0_I8 ) then
      RN_SEED0  = new_seed
    endif

    if( present(print_info) ) then
      if( print_info /= 0 ) then
        write(printseed,'(i20)') RN_SEED0
        write( iuo,1) RN_SEED0, RN_MULT, RN_ADD, RN_BITS, RN_STRIDE
        write(jtty,2) adjustl(printseed)
1       format( &
          & /,' ***************************************************', &
          & /,' * Random Number Seed       = ',i20,             ' *', &
          & /,' * Random Number Multiplier = ',i20,             ' *', &
          & /,' * Random Number Adder      = ',i20,             ' *', &
          & /,' * Random Number Bits Used  = ',i20,             ' *', &
          & /,' * Random Number Stride     = ',i20,             ' *', &
          & /,' ***************************************************',/)
2       format(' comment. using random number initial seed = ',a20)
      endif
    endif

    ! double-check on number of bits in a long int
    if( bit_size(RN_SEED)<64 ) then
      call expire( 0, 'RN_init_problem', &
        & ' ***** ERROR: <64 bits in long-int, can-t generate RN-s')
    endif
    itemp1 = 5_I8**25
    itemp2 = 5_I8**19
    itemp3 = ishft(2_I8**62-1_I8,1) + 1_I8
    itemp4 = itemp1*itemp2
    if( iand(itemp4,itemp3)/=8443747864978395601_I8 ) then
      call expire( 0, 'RN_init_problem', &
        & ' ***** ERROR: can-t do 64-bit integer ops for RN-s')
    endif

    ! set the initial particle seed
    RN_SEED  = RN_SEED0

    return
  end subroutine RN_init_problem

  !-------------------------------------------------------------------

  subroutine RN_init_particle( nps )
    ! initialize MCNP random number parameters for particle "nps"
    !
    !     * generate a new particle seed from the base seed
    !       & particle index
    !     * set the RN count to zero
    implicit none
    integer(I8), intent(in) :: nps

    RN_SEED  = RN_skip_ahead( RN_SEED0, nps*RN_STRIDE )

    return
  end subroutine RN_init_particle

  !-------------------------------------------------------------------

  subroutine expire( i, c1, c2 )
    integer,          intent(in) :: i
    character(len=*), intent(in) :: c1, c2
    write(*,*) ' ********** error: ',c1
    write(*,*) ' ********** error: ',c2
    stop i
  end subroutine expire

  !-------------------------------------------------------------------
  !###################################################################
  !#
  !#  Unit tests
  !#
  !###################################################################

  subroutine RN_test_basic
    ! test routine for basic random number generator
    implicit none
    real(R8)    :: s
    integer(I8) :: seeds(10)
    integer     :: i, j

    write(jtty,"(/,a)")  " ***** random number - basic test *****"

    ! set the seed
    call RN_init_problem( 1_I8, 1 )

    ! get the first 5 seeds, then skip a few, get 5 more - directly
    s = 0.0_R8
    do  i = 1,5
      s = s + rang()
      seeds(i) = RN_SEED
    enddo
    do  i = 6,123455
      s = s + rang()
    enddo
    do  i = 6,10
      s = s + rang()
      seeds(i) = RN_SEED
    enddo

    ! compare
    do  i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i), "  computed: ", seeds(i)
      if( seeds(i)/=RN_CHECK(i) ) then
        write(jtty,"(a)")  " ***** basic_test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_basic

  !-------------------------------------------------------------------

end module mcnp_random

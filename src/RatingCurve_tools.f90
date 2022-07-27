module RatingCurve_tools

!~**********************************************************************
!~* Purpose: Catalogue of rating curve functions
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 03/08/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. GetRCParNumber, number of parameters of the RC
!~*		2. ApplyRC, compute Q=f(H|teta)
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetRCParNumber, ApplyRC, RC_General_CheckControlMatrix,RC_General_Continuity,&
          RC_General_GetRangeFromTeta

Character(100), parameter, PUBLIC:: & ! Catalogue of available rating curves
                    RC_General="RC_General",& ! General rating curve formalism as defined by Laurent
                    RC_Linear="RC_Linear",& ! Q=a*H
                    RC_Affine="RC_Affine",& ! Q=a*H+b
                    RC_Power="RC_Power",& ! Q=a*(H-b)^c
                    RC_Kinetics_Accumulation="RC_Kinetics_Accumulation",& ! M=a*(1-exp(-kt))
                    RC_Kinetics_Elimination="RC_Kinetics_Elimination",& ! M=M0*exp(-kt)
                    RC_Kinetics_Azziz1="RC_Kinetics_Azziz1",& ! M=C*K*V*(1-exp(-kt))
                    RC_Kinetics_Azziz2="RC_Kinetics_Azziz2",& ! M=C*K*V*(1-exp(-(R/(K*V))*t))
                    RC_Kinetics_Azziz2_Mixture="RC_Kinetics_Azziz2_Mixture",& ! M=C*K1*V*(1-exp(-(R1/(K1*V))*t))+C*K2*V*(1-exp(-(R2/(K2*V))*t)) - WARNING: K1<K2 par convention
                    RC_Charbo="RC_Charbo" !a1 b1 c1 a2 b2 c2 (1=triangle,2=rectangle)

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine GetRCParNumber(RCID,npar,ControlMatrix,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the RC
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 03/08/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. RCID, ID of the rating curve
!^*		2. [ControlMatrix], optional, control Matrix in Laurent's framework
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::RCID
integer(mik), intent(in),optional::ControlMatrix(:,:)
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
integer(mik)::ncontrol

err=0;mess='';npar=undefIN

select case(RCID)
case(RC_Linear)
    npar=1
case(RC_Affine,RC_Kinetics_Accumulation,RC_Kinetics_Elimination)
    npar=2
case(RC_Power)
    npar=3
case(RC_Kinetics_Azziz1,RC_Kinetics_Azziz2)
    npar=4
case(RC_Kinetics_Azziz2_Mixture,RC_Charbo)
    npar=6
case(RC_General)
    if(.not.present(ControlMatrix)) then
        err=2;mess='GetRCParNumber: Fatal: [ControlMatrix] has to be provided for this RC';return
    endif
    ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub
    npar=3*ncontrol
case default
    err=1;mess='GetRCParNumber: Fatal: Unavailable [RCID]'
end select

end subroutine GetRCParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyRC(RCID,H,teta,ControlMatrix,Q,feas,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. RCID, ID of the rating curve
!^*		2. H, water stage
!^*		3. teta, parameters
!^*		4. [ControlMatrix], optional, control Matrix in Laurent's framework
!^* OUT
!^*		1. Q, RC-computed runoff
!^*		2. feas, feasible?
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************

character(*), intent(in)::RCID
real(mrk), intent(in)::H, teta(:)
real(mrk), intent(out)::Q
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
integer(mik), intent(in),optional::ControlMatrix(:,:)
!locals
integer(mik)::npar
real(mrk)::a1,b1,c1,a2,b2,c2,k1

err=0;mess='';feas=.true.;Q=undefIN

! Check size of teta
call GetRCParNumber(RCID=RCID,npar=npar,ControlMatrix=ControlMatrix,err=err,mess=mess)
if(err>0) then
    mess="ApplyRC: "//trim(mess);return
endif

!Apply RC
select case(RCID)
case(RC_Linear)
    Q=teta(1)*H
case(RC_Affine)
    Q=teta(1)*H+teta(2)
case(RC_Power)
    if( (H-teta(2))<0._mrk )then
        feas=.false.;return
    else
        Q=teta(1)*(H-teta(2))**teta(3)
    endif
case(RC_Kinetics_Accumulation)
    Q=teta(1) * ( 1._mrk-exp( -teta(2)*H ) )
case(RC_Kinetics_Elimination)
    Q=teta(1) * exp( -teta(2)*H ) 
case(RC_Kinetics_Azziz1)
    Q=teta(1)*teta(2)*teta(3) * ( 1._mrk-exp( -teta(4)*H ) )
case(RC_Kinetics_Azziz2)
    Q=teta(1)*teta(2)*teta(3) * ( 1._mrk-exp( -1._mrk* (teta(4)/(teta(2)*teta(3))) *H ) )
case(RC_Kinetics_Azziz2_Mixture)
    if(teta(2)>=teta(5))then
        feas=.false.;return
    else
        Q=teta(1)*teta(2)*teta(3) * ( 1._mrk-exp( -1._mrk* (teta(4)/(teta(2)*teta(3))) *H ) ) +&
          teta(1)*teta(5)*teta(3) * ( 1._mrk-exp( -1._mrk* (teta(6)/(teta(5)*teta(3))) *H ) )
    endif
case(RC_Charbo)
    if( (H-teta(2))<0._mrk)then
        feas=.false.;return
    else
        Q=teta(1)*(max(H-teta(2),0._mrk))**teta(3) & ! triangle
         -teta(1)*(max(H-teta(5),0._mrk))**teta(3) & ! triangle tronqué
         +teta(4)*(max(H-teta(5),0._mrk))**teta(6) ! rectangle
    endif
case(RC_General)
    if(.not.present(ControlMatrix)) then
        err=2;mess='ApplyRC: Fatal: [ControlMatrix] has to be provided for this RC';return
    endif
    call ApplyRC_General(H=H,teta=teta,ControlMatrix=ControlMatrix,Q=Q,feas=feas,err=err,mess=mess)
case default
    err=1;mess='ApplyRC: Fatal: Unavailable [RCID]'
end select

end subroutine ApplyRC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine RC_General_CheckControlMatrix(ControlMatrix,feas,err,mess)

!^**********************************************************************
!^* Purpose: Check Control Matrix is OK
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ControlMatrix
!^* OUT
!^*		1.feas
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::ControlMatrix(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
logical(mlk), intent(out)::feas
! locals
integer(mik)::ncontrol,nrange,i
logical(mlk)::mask(size(ControlMatrix,dim=2))
err=0;mess='';feas=.true.

ncontrol=size(ControlMatrix,dim=2)
nrange=size(ControlMatrix,dim=1)

if(ncontrol/=nrange) then
    err=1;mess='RC_General_CheckControlMatrix:FATAL:ncontrol/=nrange';feas=.false.;return
endif

if (.not. all( ControlMatrix==0 .or. ControlMatrix==1) ) then
    err=2;mess='RC_General_CheckControlMatrix:FATAL:Some elements of ControlMatrix differ from 0/1';feas=.false.;return
endif

if(nrange>1) then ! check matrix is lower-diagonal
    do i=1,nrange-1
        if( any(ControlMatrix(i,i+1:)>0) ) then
            err=3;mess='RC_General_CheckControlMatrix:FATAL:ControlMatrix should be lower-diagonal';feas=.false.;return
        endif
    enddo
endif

if(nrange>1) then ! check each change of range activates exactly one control
    do i=2,nrange
        mask=ControlMatrix(i,:) - ControlMatrix(i-1,:)>0
        if (count(mask)/=1) then
            err=4;mess='RC_General_CheckControlMatrix:FATAL: Exactly one control should be activated when changing range';feas=.false.;return
        endif
    enddo
endif

do i=1,nrange ! check that change to the kth range activates the kth control
    if(ControlMatrix(i,i)/=1) then
        err=5;mess='RC_General_CheckControlMatrix:FATAL: change to the kth range should activate the kth control';feas=.false.;return
    endif
enddo

end subroutine RC_General_CheckControlMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine RC_General_Continuity(a,c,k,ControlMatrix,b,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply Continuity constraints & evaluate feasability
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.a
!^*		2.c
!^*		3.k
!^*		4.ControlMatrix
!^* OUT
!^*		1.feas, feasability
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.b
!^**********************************************************************

real(mrk), intent(in)::a(:),c(:),k(:)
integer(mik),intent(in):: ControlMatrix(:,:)
real(mrk), intent(inout)::b(:)
logical(mlk), intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::ncontrol,i,j
real(mrk)::som

err=0;mess='';feas=.true.

ncontrol=size(a)

if(a(1)<=0._mrk .or. c(1)==0._mrk) then 
    feas=.false.;return
endif
if(ncontrol<=1) return ! nothing to do - b already specified in the single-range case

do j=2,ncontrol
    if(j<ncontrol) then
        if(k(j)<=k(j-1)) then
            feas=.false.;return
        endif
    endif
    ! feasability checks
    if( any( k(j-1)<b(1:j-1) ) .or. a(j)<=0._mrk .or. c(j)==0._mrk) then
        feas=.false.;return
    endif
    som=0._mrk
    do i=1,j-1
        som=som+real(ControlMatrix(j-1,i)-ControlMatrix(j,i),mrk)*a(i)*(k(j-1)-b(i))**c(i)
    enddo
    if(som<0._mrk) then
        feas=.false.;return
    endif
    b(j)=k(j-1)-((1._mrk/a(j))*som)**(1._mrk/c(j))
enddo

end subroutine RC_General_Continuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RC_General_GetRangeFromTeta(H,teta,ControlMatrix,range,err,mess)

!^**********************************************************************
!^* Purpose: find the range where H belongs, given theta & the control matrix
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:25/04/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. H
!^*		2. teta
!^*		2. ControlMatrix
!^* OUT
!^*		1.range
!^*		2.err
!^*		3.mess
!^**********************************************************************

real(mrk), intent(in)::H,teta(:)
integer(mik), intent(in)::ControlMatrix(:,:)
integer(mik),intent(out)::range,err
character(*), intent(out)::mess
!locals
integer(mik)::ncontrol,i,nk,p
real(mrk), allocatable::k(:)

err=0;mess='';range=undefIN;

! Get number of controls (=number of ranges)
ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub

! check size of teta
if(size(teta)/=3*ncontrol) then
    err=1;mess='RC_General_GetRangeFromTheta: Fatal: Incorrect size [theta]';return
endif

! Allocate and assign values to k
if(ncontrol<=1) then
    range=1;return
else
    if(allocated(k)) deallocate(k);allocate(k(ncontrol-1))
    do i=2,ncontrol
        p=3*(i-1)
        k(i-1)=teta(p+1);
    enddo
endif

nk=ncontrol-1
do i=1,nk
    If(H<k(i)) then
        range=i;return
    else
        cycle
    endif      
enddo

range=nk+1

end subroutine RC_General_GetRangeFromTeta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private subs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyRC_General(H,teta,ControlMatrix,Q,feas,err,mess)

!^**********************************************************************
!^* Purpose: General derivation of RC based on controls and ranges 
!^* (cf Laurent's formalization)
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. H, stage
!^*		2. teta, RC parameters
!^*		3. ControlMatrix, control Matrix in Laurent's framework
!^* OUT
!^*		1.Q, runoff
!^*		2.feas, feasability
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************

real(mrk), intent(in)::H, teta(:)
real(mrk), intent(out)::Q
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
integer(mik), intent(in)::ControlMatrix(:,:)
!locals
integer(mik)::ncontrol,i,p,range
real(mrk), allocatable::a(:),b(:),c(:),k(:)

err=0;mess='';Q=undefRN;feas=.true.

! Get number of controls (=number of ranges)
ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub

! check size of teta
if(size(teta)/=3*ncontrol) then
    err=1;mess='ApplyRC_General: Fatal: Incorrect size [teta]';feas=.false.;return
endif

! Allocate and assign values to a,b,c,k
allocate(a(ncontrol),b(ncontrol),c(ncontrol),k(ncontrol-1))
a(1)=teta(1);b(1)=teta(2);c(1)=teta(3)
if(ncontrol>1) then
    do i=2,ncontrol
        p=3*(i-1)
        k(i-1)=teta(p+1);a(i)=teta(p+2);c(i)=teta(p+3)
    enddo
endif

! Apply Continuity Constraint & check feasability
call RC_General_Continuity(a,c,k,ControlMatrix,b,feas,err,mess)
if(.not.feas) return

! Retrieve range corresponding to stage H
range=RC_General_GetRange(H,k)

! Apply equation
Q=0._mrk
do i=1,range
    if(ControlMatrix(range,i)==1) then
        if((H-b(i))<0._mrk) then
            Q=0._mrk;return
        endif
        Q=Q+a(i)*(H-b(i))**c(i)
    else
        if((H-b(i))<0._mrk) then
            write(*,*) 'WARNING: You should NEVER EVER arrive there!!!'
            write(*,*) 'If you see this message, please contact developpers !!!'
            feas=.false.;return
        endif
    endif
enddo

end subroutine ApplyRC_General
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function RC_General_GetRange(H,k)

!#**********************************************************************
!#* Purpose: find the range where H belongs
!#**********************************************************************
!#* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*		1. H
!#*		2. k
!#* OUT
!#*		1.RC_General_GetRange
!#**********************************************************************

real(mrk), intent(in)::H,k(:)
integer(mik)::RC_General_GetRange
!locals
integer(mik)::nk,i

nk=size(k)
do i=1,nk
    If(H<k(i)) then
        RC_General_GetRange=i;return
    else
        cycle
    endif      
enddo

RC_General_GetRange=nk+1

end function RC_General_GetRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module RatingCurve_tools
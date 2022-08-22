module BaRatin_tools

!~**********************************************************************
!~* Purpose: compute prior/post for Bayesian Rating curve fitting (BaRatin)
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified:26/09/2013
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1.
!~*     2.
!~*     3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use BayesianEstimation_tools, only:PriorListType

implicit none
Private
public :: &! main subroutines to handle the probabilistic model behind BaRatin
          BaRatin_Fit, BaRatin_PostProcess,Sigmafunk_GetParNumber,Sigmafunk_Apply,&
          BaRatin_Propagation,BaRatin_LoadRCobject,BaRatin_H2Q,&
          ! Utilities to read configuration files
          BaRatinConfig_Read_Workspace,&
          BaRatinConfig_Read_Data,BaRatinConfig_Read_MCMC,&
          BaRatinConfig_Read_RatingCurve,BaRatinConfig_Read_ControlMatrix,&
          BaRatinConfig_Read_RemnantSigma,BaRatinConfig_Read_RunOptions,&
          BaRatinConfig_Read_PriorRC,BaRatinConfig_Read_PostProcessing,&
          BaRatinConfig_Read_H2QPropagation,BaRatinConfig_Read_Aggregation,&
          ! Reading other files
          BaRatin_ReadData,&
          ! Messages & errors utilities
          BaRatin_ConsoleMessage,BaRatin_Fatal_Exit

! variables globally available to this module
type, public:: RCType ! the "rating curve" object
    character(100)::RCID='AintGotNoName' ! RC formula - see RatingCurve_tools for the catalogue
    integer(mik)::nteta=undefIN ! number of parameters for the rating curve
    integer(mik)::nobs=undefIN ! number of observations used to fit the RC
    logical:: Herror=.false. ! is any stage value affected by an error? if not, this is WLS-like fitting
    integer(mik)::nHerror=undefIN ! number of H values affected by an error
    real(mrk), allocatable::Hobs(:),Qobs(:) ! observed series of H/Q
    real(mrk), allocatable::Hsigma(:),Qsigma(:) ! series of H/Q standard errors
    type(PriorListType), allocatable:: PriorList_RemnantSigma(:) ! priors for remnant std - see BayesianEstimation_tools module for PriorListType definition
    type(PriorListType), allocatable:: PriorList_teta(:) ! priors for RC parameters
    integer(mik),allocatable::ControlMatrix(:,:)
    integer(mik)::ncontrol=1 ! [only for General-RC] number of hydraulic controls
    integer(mik)::nsigmaf=1 ! number of parameters for remnant sdev.
    character(250)::RemnantSigma_funk=undefCH
    real(mrk)::mv=-9999._mrk
end type RCType
type(RCType),public:: RC

integer(mik),parameter::RemnantOption_def=1 ! 0=systematic,1=random
integer(mik),parameter::hcol=1,tcol=2,rcol=3 ! meaning of columns for propagation matrix
integer(mik),parameter,public::messID_Read=-3,messID_Open=-2,messID_Write=-8
character(5),parameter,public::BaRatin_realFmt='e15.6' ! 2DO: scan this module and replace all hard-coded occurrences

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_LoadRCobject(Hobs,Hsigma,Qobs,Qsigma,& ! observations and uncertainties
                        RCID,&               ! Rating Curve ID
                        ControlMatrix,& ! optional Control Matrix for Laurent's general formalisation
                        RemnantSigma_funk,& ! chosen f in sdev=f(Qrc)
                        PriorList_teta,PriorList_RemnantSigma,&
                        err,mess)! error handling
!^**********************************************************************
!^* Purpose: Load RC object
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. Hobs, observed stage
!^*     2. Hsigma, std of stage error
!^*     3. Qobs, observed runoff
!^*     4. Qsigma, std of runoff error
!^*     5. RCID, ID of the rating curve
!^*     6. [ControlMatrix], optional Control Matrix for Laurent's general formalisation
!^*     7. RemnantSigma_funk, function for remnant sigma
!^*     8. PriorList_teta, priors for RC parameters
!^*     9. PriorList_RemnantSigma, priors for remnant sigma parameters
!^* OUT
!^*     1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2.mess, error message
!^**********************************************************************
use RatingCurve_tools, only:GetRCParNumber
real(mrk), intent(in)::Hobs(:),Hsigma(:),Qobs(:),Qsigma(:)
character(*), intent(in)::RCID
integer(mik),intent(in),optional:: ControlMatrix(:,:)
Type(PriorListType), intent(in)::PriorList_teta(:), PriorList_RemnantSigma(:)
character(*),intent(in)::RemnantSigma_funk
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatin_LoadRCobject'
integer(mik)::nHerror

err=0;mess=''
! Size checks
if( size(Hobs)/=size(Hsigma) .or. size(Hobs)/=size(Qobs) .or. size(Hobs)/=size(Qsigma) .or. &
    size(Hsigma)/=size(Qobs) .or. size(Hsigma)/=size(Qsigma) .or. size(Qobs)/=size(Qsigma) ) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(20));return
endif
! Populate object
RC%RCID=RCID
call GetRCParNumber(RCID=RCID,ControlMatrix=ControlMatrix,npar=RC%nteta,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
RC%nobs=size(Hobs)
if(allocated(RC%Hobs)) deallocate(RC%Hobs);allocate(RC%Hobs(RC%nobs))
RC%Hobs=Hobs
if(allocated(RC%Hsigma)) deallocate(RC%Hsigma);allocate(RC%Hsigma(RC%nobs))
RC%Hsigma=Hsigma
if(allocated(RC%Qsigma)) deallocate(RC%Qsigma);allocate(RC%Qsigma(RC%nobs))
RC%Qsigma=Qsigma
if(allocated(RC%Qobs)) deallocate(RC%Qobs);allocate(RC%Qobs(RC%nobs))
RC%Qobs=Qobs
RC%Herror = any(Hsigma/=0._mrk)
nHerror=count(RC%Hsigma /= 0._mrk)
RC%nHerror=nHerror
! priors
if(allocated(RC%PriorList_teta)) deallocate(RC%PriorList_teta);allocate(RC%PriorList_teta(RC%nteta))
RC%PriorList_teta=PriorList_teta
if(allocated(RC%PriorList_RemnantSigma)) then
    deallocate(RC%PriorList_RemnantSigma)
    allocate(RC%PriorList_RemnantSigma(size(PriorList_RemnantSigma)))
endif
RC%PriorList_RemnantSigma=PriorList_RemnantSigma
! control matrix
if(present(ControlMatrix)) then
    if(allocated(RC%ControlMatrix)) deallocate(RC%ControlMatrix)
    allocate(RC%ControlMatrix(size(ControlMatrix,dim=1),size(ControlMatrix,dim=2)))
    RC%ControlMatrix=ControlMatrix
    RC%ncontrol=size(ControlMatrix,dim=1)
else
    RC%ncontrol=1
endif
! handle RemnantSigma
if(present(ControlMatrix)) then
    RC%ncontrol=size(ControlMatrix,dim=1)
else
    RC%ncontrol=undefIN
endif
RC%RemnantSigma_funk=RemnantSigma_funk
call Sigmafunk_GetParNumber(RC%RemnantSigma_funk, RC%nsigmaf, err, mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine BaRatin_LoadRCobject

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetLogPost(teta,RemnantSigma,Htrue,& ! inferred quantities
                      RC,& ! rating curve object
                      lp, feas, isnull,err,mess) ! outputs

!^**********************************************************************
!^* Purpose: compute post(teta,RemnantSigma,Htrue)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 03/10/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. teta, RC parameters
!^*     2. RemnantSigma, std of remnant errors
!^*     3. Htrue, estimated "true" stage
!^*     4. RC, rating curve object
!^* OUT
!^*     1.lp, log-posterior
!^*     2.feas, feasible?
!^*     3.is null, is (natural) posterior = zero?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************
use BayesianEstimation_tools, only:GetLogPrior, PriorListType
use RatingCurve_tools, only:GetRCParNumber, ApplyRC,RC_General_GetRangeFromTeta
use Distribution_tools, only: GetPdf, GAUSS

real(mrk), intent(in)::teta(:),RemnantSigma(:),Htrue(:)
type(RCType),intent(in):: RC
real(mrk), intent(out)::lp
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetLogPost'
real(mrk)::prior_RemnantSigma, prior_teta, prior, Qhat, Qlik, Hlik, logp, mu, sig
integer(mik)::n, npar, i, ncontrol,range,j,SigmaFunk_np

!Init
err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN
n=size(RC%Hobs)

! Size checks
if( size(Htrue)/=size(RC%Hobs)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(15));return
endif

if(size(RC%PriorList_RemnantSigma)/=size(RemnantSigma)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(16));return
endif

if(size(RC%PriorList_teta)/=size(teta)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(17));return
endif

call GetRCParNumber(RCID=RC%RCID,npar=npar,ControlMatrix=RC%ControlMatrix,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(npar/=size(teta)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(18));return
endif

call Sigmafunk_GetParNumber(RC%RemnantSigma_funk, SigmaFunk_np, err, mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if(size(RemnantSigma)/=SigmaFunk_np) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(19));return
endif

! Get Prior
call GetLogPrior(teta=RemnantSigma,PriorList=RC%PriorList_RemnantSigma,&
                lp=prior_RemnantSigma, feas=feas, isnull=isnull,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return

call GetLogPrior(teta=teta,PriorList=RC%PriorList_teta,&
                lp=prior_teta, feas=feas, isnull=isnull,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return
prior=prior_teta + prior_RemnantSigma

! Get Q-likelihood
Qlik=0._mrk
do i=1,n
    call ApplyRC(RCID=RC%RCID,H=Htrue(i),teta=teta,ControlMatrix=RC%ControlMatrix,Q=Qhat,feas=feas,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(.not. feas) return
    call Sigmafunk_Apply(RC%RemnantSigma_funk, RemnantSigma, Qhat, sig, err, mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
    if(sig<=0._mrk) then;
        feas=.false.;return
    else
        sig=sqrt(RC%Qsigma(i)**2+sig**2)
    endif
    ! Qobs ~ N(mu = Qhat, sigma² = sigma_Q² + sigma_remnant²)
    mu=Qhat
    call GetPdf(DistId=GAUSS,x=RC%Qobs(i),par=(/mu,sig/),loga=.true.,pdf=logp,&
                feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull) return
    Qlik=Qlik+logp
enddo

! Get H-likelihood
Hlik=0._mrk
do i=1,n
    if(RC%Hsigma(i) == 0._mrk) then ! no stage error
        if(Htrue(i) == RC%Hobs(i)) then !OK
            cycle
        else ! not OK...
            isnull=.true.;return
        endif
    endif
    call GetPdf(DistId=GAUSS,x=RC%Hobs(i),par=(/Htrue(i),RC%Hsigma(i)/),loga=.true.,pdf=logp,&
                feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull) return
    Hlik=Hlik+logp
enddo

! Get Posterior
lp = prior + Hlik + Qlik

end subroutine GetLogPost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Posterior_wrapper(x,feas,isnull,fx,fAux,err,mess)

!^**********************************************************************
!^* Purpose: wrapper to GetLogPost to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:04/08/2010
!^**********************************************************************
!^* Comments: mv handled out of here
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.x, teta
!^* OUT
!^*     1. see GetLogPost
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^* INOUT
!^*     1.
!^*     2.
!^*     3.
!^**********************************************************************

real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Posterior_wrapper'
real(mrk), allocatable:: Htrue(:)
integer(mrk)::k, nHerror, i

if( .not. RC%Herror ) then
    call GetLogPost(teta=x(1:RC%nteta),RemnantSigma=x( (RC%nteta+1):(RC%nteta+RC%nsigmaf) ),Htrue=RC%Hobs,&
                    RC=RC,&
                    lp=fx, feas=feas, isnull=isnull,err=err,mess=mess)
else
    nHerror=count(RC%Hsigma /= 0._mrk)
    if( size(x)/=(RC%nteta+RC%nsigmaf+nHerror) ) then
        err=1;mess=trim(procname)//':'//trim(BaRatin_message(14));return
    endif
    if(allocated(Htrue)) deallocate(Htrue);allocate(Htrue(RC%nobs))
    Htrue=RC%Hobs;k=0
    do i=1,RC%nobs
        if(RC%Hsigma(i) /= 0._mrk) then
            k=k+1;Htrue(i)=x(RC%nteta+RC%nsigmaf+k)
        endif
    enddo

    call GetLogPost(teta=x(1:RC%nteta),RemnantSigma=x( (RC%nteta+1):(RC%nteta+RC%nsigmaf) ),Htrue=Htrue,&
                    RC=RC,&
                    lp=fx, feas=feas, isnull=isnull,err=err,mess=mess)
endif

if(present(fAux)) fAux=UndefRN

end subroutine Posterior_wrapper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_Fit(RC,& ! Rating curve object
                       !!!!!!! Tuning of the MCMC sampler !!!!!!!!
                       teta0,RemnantSigma0, &!initial values for teta and remnant std
                       teta_std0,RemnantSigma_std0,& ! initial values for the std of jump distribution
                       nAdapt,nCycles,&
                       MinMoveRate,MaxMoveRate,&
                       DownMult,UpMult,&
                       !!!!!!! END Tuning of the MCMC sampler !!!!!!!!
                       OutFile, & ! Output file (for MCMC samples)
                       err,mess)! error handling


!^**********************************************************************
!^* Purpose: Performs MCMC sampling
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
!^*     1.RC, rating curve object
!^*     2 and above: properties of the MCMC sampler
!^*     3. OutFile, Output file (for MCMC samples)
!^* OUT
!^*     1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2.mess, error message
!^**********************************************************************
use MCMCStrategy_tools, only:Adaptive_Metro_OAAT
use utilities_dmsl_kit, only:number_string

! INPUTS
type(RCType),intent(in):: RC
real(mrk), intent(in)::teta0(:),RemnantSigma0(:)
real(mrk), intent(in)::teta_std0(:),RemnantSigma_std0(:)
integer(mik), intent(in)::nAdapt,nCycles
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in)::OutFile
! OUTPUTS
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaRatin_Fit'
real(mrk)::fx
integer(mik)::nHerror, i, k
real(mrk), allocatable:: start(:), startStd(:)
character(100), allocatable:: headers(:)

err=0;mess=''

! Size checks
if( size(teta0)/=size(teta_std0) .or. size(teta0)/=size(RC%PriorList_teta) .or. &
    size(teta_std0)/=size(RC%PriorList_teta)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_message(12));return
endif
if( size(RemnantSigma0)/=size(RemnantSigma_std0) .or. size(RemnantSigma0)/=size(RC%PriorList_RemnantSigma) .or. &
    size(RemnantSigma_std0)/=size(RC%PriorList_RemnantSigma)) then
    err=1;mess=trim(procname)//':'//trim(Baratin_message(13));return
endif

! MCMC sampling
nHerror=count(RC%Hsigma /= 0._mrk)
if(allocated(start)) deallocate(start)
allocate(start(RC%nteta+RC%nsigmaf+nHerror))
if(allocated(startStd)) deallocate(startStd)
allocate(startStd(RC%nteta+RC%nsigmaf+nHerror))
if(allocated(headers)) deallocate(headers)
allocate(headers(RC%nteta+RC%nsigmaf+nHerror+1))
start(1:(RC%nteta+RC%nsigmaf))=(/teta0,RemnantSigma0/)
startStd(1:(RC%nteta+RC%nsigmaf))=(/teta_std0,RemnantSigma_std0/)
do i=1,RC%nteta
    headers(i)="Teta"//trim(number_string(i))
enddo
do i=1,RC%nsigmaf
    headers(RC%nteta+i)="RemnantSTD"//trim(number_string(i))
enddo
if(RC%Herror) then
    k=0
    do i=1,RC%nobs
        if(RC%Hsigma(i) /= 0._mrk) then
            k=k+1
            start(RC%nteta+RC%nsigmaf+k)=RC%Hobs(i)
            startStd(RC%nteta+RC%nsigmaf+k)=RC%Hsigma(i)
            headers(RC%nteta+RC%nsigmaf+k)="trueH"//trim(number_string(i))
        endif
    enddo
endif
headers(RC%nteta+RC%nsigmaf+nHerror+1)="LogPost"

! GO!
call Adaptive_Metro_OAAT(f=Posterior_wrapper,x=start,&
                fx=fx,std=startStd,&
                nAdapt=nAdapt,nCycles=nCycles,&
                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                DownMult=DownMult,UpMult=UpMult,&
                OutFile=OutFile,headers=headers, err=err,mess=mess)

end subroutine BaRatin_Fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_ReadMCMC(MCMCFile,Nburn,Nread,Nslim,& ! Read properties
                        mcmc,err,mess)! error handling
!^**********************************************************************
!^* Purpose: Load MCMC file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. MCMCFile, containing MCMC samples
!^*     2. Nburn, number of lines to discard
!^*     3. Nread, number of lines to read
!^*     4. Nslim, only one line every Nslim will be used
!^* OUT
!^*        1.mcmc
!^*        1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use RatingCurve_tools, only:RC_General,RC_General_Continuity
character(*), intent(in)::MCMCFile
integer(mik),intent(in)::Nburn,Nread,Nslim
integer(mik), intent(out)::err
character(*),intent(out)::mess
real(mrk),pointer::mcmc(:,:)
! locals
character(250),parameter::procname='BaRatin_ReadMCMC'
integer(mik)::unt,i,j,k,n,errcode,nMCMC,p
real(mrk), allocatable::Temp(:,:)
real(mrk), allocatable::aRC(:,:),bRC(:,:),cRC(:,:),kRC(:,:)
logical::feas

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

open(unit=unt,file=trim(MCMCFile), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(MCMCFile));endif

read(unt,*,iostat=err) !headers
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(MCMCFile));endif
do i=1,Nburn ! burn
    read(unt,*,iostat=errcode)
    if(errcode/=0) then ! problem reading file
        err=1;mess=trim(procname)//':'//trim(BaRatin_Message(2));return
    endif
enddo
if(allocated(Temp)) deallocate(Temp);allocate(Temp(Nread,RC%nteta+RC%nsigmaf+RC%nHerror+1))
do i=1,Nread ! read used lines
    k=i
    read(unt,*,iostat=errcode) Temp(i,:)
    if(errcode/=0) then ! problem reading file
        err=-1;k=k-1
        mess=trim(procname)//':'//trim(BaRatin_Message(3));exit
    endif
enddo
close(unt)
! Slim
n=size(Temp(1:k:nSlim,:),dim=1) ! k/Nslim

! Extend mcmc table for RC of type RC_General
if(RC%RCID==RC_General) then
    if(associated(mcmc)) nullify(mcmc);allocate(mcmc(n,RC%nteta+RC%nsigmaf+RC%nHerror+1+RC%Ncontrol-1)) ! Add computed b's to all result files
    mcmc(:,1:(RC%nteta+RC%nsigmaf+RC%nHerror+1))=Temp(1:k:nSlim,:)
    nMCMC=size(mcmc,dim=1)
    ! retrieve a's, k's and c's
    allocate(aRC(nMCMC,RC%ncontrol),bRC(nMCMC,RC%ncontrol),cRC(nMCMC,RC%ncontrol),kRC(nMCMC,RC%ncontrol-1))
    do j=1,nMCMC
        aRC(j,1)=mcmc(j,1);bRC(j,1)=mcmc(j,2);cRC(j,1)=mcmc(j,3)
        if(RC%ncontrol>1) then
            do i=2,RC%ncontrol
                p=3*(i-1)
                kRC(j,i-1)=mcmc(j,p+1);aRC(j,i)=mcmc(j,p+2);cRC(j,i)=mcmc(j,p+3)
            enddo
        endif
        call RC_General_Continuity(a=aRC(j,:),c=cRC(j,:),k=kRC(j,:),ControlMatrix=RC%ControlMatrix,b=bRC(j,:),&
                                   feas=feas,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    enddo
    mcmc(:,(RC%nteta+RC%nsigmaf+RC%nHerror+2):(RC%nteta+RC%nsigmaf+RC%nHerror+1+RC%ncontrol-1))=bRC(:,2:)
else
    if(associated(mcmc)) nullify(mcmc);allocate(mcmc(n,RC%nteta+RC%nsigmaf+RC%nHerror+1))
    mcmc=Temp(1:k:nSlim,:)
endif

end subroutine BaRatin_ReadMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_H2Q(RC,& ! Rating curve object
                       Hx,Hx_Std,Hx_Bias,Hx_BiasIndx,&
                       mcmc,nsim,PropagateWhat,&
                       PrintCounter,RemnantOption,&
                       SpagFilePrefix,EnvelopFilePrefix,&
                       SaveSpag,SaveEnvelop,&
                       err,mess)! error handling
!^**********************************************************************
!^* Purpose: Generates Q spaghettis
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        0. RC, rating curve object
!^*        1. Hx, H values
!^*        2. [Hx_Std], optional: uncertainty in Hs
!^*        3. [Hx_Bias], optional: standard deviation of the H-bias term
!^*        4. [Hx_BiasIndx], optional: index for bias-resampling
!^*        5. [mcmc], mcmc sample. If not present, will generate teta from the prior
!^*        6. [nsim], only used for prior sampling, not used is mcmc is provided
!^*        7. PropagateWhat, 0-1 vector of size 3, propagate h-U, teta-U, remnant-U?
!^*        8. [PrintCounter], print counter in console during iterations? default .false.
!^*        9. [RemnantOption], 0=systematic,1=independent,default: see RemnantOption_def
!^*        10. SpagFilePrefix, prefix for spaghetti file.
!^*        11. EnvelopFilePrefix, prefix for envelop file.
!^*        12. SaveSpag, save spag file?
!^*        13. SaveEnvelop, save envelop file?
!^* OUT
!^*        1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        2.mess, error message
!^**********************************************************************
use RatingCurve_tools, only:ApplyRC
use Distribution_tools, only:Generate,GenerateSample,GAUSS,GetMode
use utilities_dmsl_kit,only:getSpareUnit,number_string

type(RCType),intent(in):: RC
real(mrk), intent(in)::Hx(:)
real(mrk), intent(in), optional::Hx_Std(:),Hx_Bias(:),mcmc(:,:)
integer(mik), intent(in), optional::Hx_BiasIndx(:)
integer(mik),intent(in)::PropagateWhat(3)
integer(mik),intent(in),optional::nsim
integer(mik),intent(in),optional::RemnantOption
logical, intent(in), optional::PrintCounter
logical,intent(in)::SaveSpag,SaveEnvelop
character(*), intent(in)::SpagFilePrefix,EnvelopFilePrefix
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250), dimension(9),parameter::head_envelop=(/"H              ",&
                                                       "Q_Modal        ",&
                                                       "Q_Median       ",&
                                                       "Q_q2.5         ",&
                                                       "Q_q97.5        ",&
                                                       "Q_q16          ",&
                                                       "Q_q84          ",&
                                                       "Q_MaxPost-sig  ",&
                                                       "Q_MaxPost+sig  "/)
character(250),parameter::procname='BaRatin_H2Q'
integer(mik),parameter::Nrefresh=1
integer(mik)::NHx,Nmcmc,i,j,k,Ropt,unt,unt2,ml(1),next,Nb
real(mrk)::dev,res,h,teta(RC%nteta),gamma(RC%nsigmaf),Mteta(RC%nteta),Mgamma(RC%nsigmaf),env(8)
logical::feas,printC
character(250),allocatable::head(:)
real(mrk),allocatable::eps(:),spag(:),priorSample(:,:),b(:,:)
character(250):: fmt

err=0;mess=''
! basic checks
if(PropagateWhat(rcol)==1 .and. (.not.present(mcmc))) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(7));return
endif

if( (.not.present(Hx_Std)).and.PropagateWhat(hcol)==1) then
    err=1;mess=trim(procname)//':'//trim(Baratin_Message(8));return
endif
if(present(mcmc)) then
    Nmcmc=size(mcmc,dim=1)
else
    if(.not.present(nsim)) then
        err=1;mess=trim(procname)//':'//trim(Baratin_Message(9));return
    endif
    Nmcmc=Nsim
endif
if(allocated(eps)) deallocate(eps);allocate(eps(Nmcmc))
if(allocated(spag)) deallocate(spag);allocate(spag(Nmcmc))
if(allocated(head)) deallocate(head);allocate(head(Nmcmc+1))

! handle defaults
if(present(RemnantOption)) then;Ropt=RemnantOption;else;Ropt=RemnantOption_def;endif
if(present(PrintCounter)) then;printC=PrintCounter;else;printC=.false.;endif
NHx=size(Hx);

! prepare result files
if(SaveSpag) then
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unt,file=trim(SpagFilePrefix)//'_'//&
                &trim(number_string(PropagateWhat(hcol)))//&
                &trim(number_string(PropagateWhat(tcol)))//&
                &trim(number_string(PropagateWhat(rcol)))//'.txt',status='replace')
    head(1)='H'
    do i=1,Nmcmc;head(i+1)="Q_"//trim(number_string(i));enddo
    fmt='('//trim(number_string(Nmcmc+1))//'A15)'
    write(unt,trim(fmt)) head
endif
if(SaveEnvelop) then
    call getSpareUnit(unt2,err,mess)
    if(err/=0) then;mess='BaRatin_H2Q:'//trim(mess);return;endif
    open(unt2,file=trim(EnvelopFilePrefix)//'_'//&
                &trim(number_string(PropagateWhat(hcol)))//&
                &trim(number_string(PropagateWhat(tcol)))//&
                &trim(number_string(PropagateWhat(rcol)))//'.txt',status='replace')
    write(unt2,'(9A15)') head_envelop
endif

! Go!!

! Generate Remnant errors here if remnant error is systematic
if(Ropt==0 .and. PropagateWhat(rcol)==1) then
    call GenerateSample(DistId=GAUSS,par=(/0._mrk,1._mrk/),&
                        gen=eps,feas=feas,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
endif

! Get modal parameters
if(present(mcmc)) then ! posterior mode
    ml=maxloc(mcmc(:,RC%nteta+RC%nsigmaf+RC%nHerror+1))
    Mteta=mcmc(ml(1),1:RC%nteta)
    Mgamma=mcmc(ml(1),(RC%nteta+1):(RC%nteta+RC%nsigmaf))
else ! prior mode
    do k=1,RC%nteta
        call GetMode(DistId=RC%PriorList_teta(k)%dist,par=RC%PriorList_teta(k)%par,&
                 m=Mteta(k),feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) then
            err=1;mess=trim(procname)//':'//trim(Baratin_Message(10));return
        endif
    enddo
    Mgamma=undefRN
endif

! generate prior sample if required
if(.not.present(mcmc)) then
    if(allocated(priorSample)) deallocate(priorSample);allocate(priorSample(Nmcmc,RC%nteta))
    do k=1,RC%nteta
        call generateSample(DistId=RC%PriorList_teta(k)%dist,par=RC%PriorList_teta(k)%par,&
                            gen=priorSample(:,k),feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) then
            err=1;mess=trim(procname)//':'//trim(Baratin_Message(11));return
        endif
    enddo
endif

! prepare normalized systematic bias values
if(PropagateWhat(hcol)==1) then
    Nb=maxval(Hx_BiasIndx)
    if(allocated(b)) deallocate(b);allocate(b(Nmcmc,Nb))
    do i=1,Nmcmc
        call GenerateSample(DistId=GAUSS,par=(/0._mrk,1._mrk/),&
                      gen=b(i,:),feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    enddo
endif
! Start big loop
do j=1,NHx ! loop on stage values
    if(printC) then
        if(j==1) next=Nrefresh
        if(j>=0.01*next*Nhx) then
!        if(mod( 100.*real(j)/real(NHx), 1.)==0.) then
            write(*,'(I4,A)') next,'% DONE'
            next=next+Nrefresh
        endif
    endif
    do i=1,Nmcmc ! loop on teta values
        if(PropagateWhat(hcol)==1) then
            if(Hx_Std(j)>0._mrk) then
                call Generate(DistId=GAUSS,par=(/Hx(j),Hx_Std(j)/),&
                          gen=h,feas=feas,err=err,mess=mess)
                if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            else
                h=Hx(j)
            endif
            if(Hx_BiasIndx(j)>0) then
                h=h+b(i,Hx_BiasIndx(j))*Hx_Bias(j)
            endif
        else
            h=Hx(j)
        endif
        if(PropagateWhat(tcol)==0) then ! use modal teta
            teta=Mteta
            gamma=Mgamma
        else ! propagate teta uncertainty
            if(present(mcmc)) then ! use mcmc sample from posterior
                teta=mcmc(i,1:RC%nteta)
                gamma=mcmc(i,(RC%nteta+1):(RC%nteta+RC%nsigmaf))
            else ! sample from prior
                teta=priorSample(i,1:RC%nteta)
                gamma=undefRN
            endif
        endif
        ! apply RC
        Call ApplyRC(RCID=RC%RCID,H=h,teta=teta,&
                    ControlMatrix=RC%ControlMatrix,Q=spag(i),&
                    feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) then;spag(i)=RC%mv;cycle;endif
        ! add remnant error if required
        if(PropagateWhat(rcol)==1) then
            call Sigmafunk_Apply(funk=RC%RemnantSigma_funk, par=gamma, &
                            Qrc=spag(i), res=res, err=err, mess=mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
            if(Ropt==1) then ! generate a new remnant error for every new h
                call Generate(DistId=GAUSS,par=(/0._mrk,1._mrk/),&
                          gen=dev,feas=feas,err=err,mess=mess)
                if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            else
                dev=eps(i) ! use the same normalized error for all h's
            endif
        else
            dev=0._mrk;res=0._mrk
        endif
        spag(i)=spag(i)+dev*res
    enddo
    if(SaveSpag) then
        ! write spaghetti to file
        fmt='('//trim(number_string(Nmcmc+1))//'e15.6)'
        write(unt,trim(fmt)) Hx(j),spag
    endif
    if(SaveEnvelop) then
        ! Get envelop and write it
        call BaRatin_GetEnvelop(RC,Hx(j),Mteta,spag,env,err,mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        write(unt2,'(9e15.6)') Hx(j),env
    endif
enddo

if(SaveSpag) close(unt)
if(SaveEnvelop) close(unt2)
deallocate(spag);deallocate(eps)

end subroutine BaRatin_H2Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_GetEnvelop(RC,Hx,mode,Spag,Env,err,mess)
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile
use RatingCurve_tools, only:ApplyRC
use numerix_dmsl_kit, only:quicksort

real(mrk), intent(in)::Hx,mode(:),Spag(:)
type(RCType), intent(in)::RC
real(mrk),intent(out)::Env(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatin_GetEnvelop'
integer(mik)::i
logical::feas
real(mrk)::std,arr(size(Spag))

! Max-Post RC
Call ApplyRC(RCID=RC%RCID,H=Hx,teta=mode,&
                ControlMatrix=RC%ControlMatrix,Q=Env(1),&
                feas=feas,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.feas) Env(1)=RC%mv
! Envelops
! Sort outside of quantiles subs to decrease computing time
arr=Spag;call quicksort(arr=arr,ascnd=.true.,err=err)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetEmpiricalStats(x=Spag,std=std,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(any(Spag==RC%mv)) then
    Env(7)=RC%mv;Env(8)=RC%mv
else
    Env(7)=Env(1)-std;Env(8)=Env(1)+std
endif
call GetEmpiricalQuantile(p=0.5_mrk,x=arr,IsXSorted=.true.,q=Env(2),err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetEmpiricalQuantile(p=0.16_mrk,x=arr,IsXSorted=.true.,q=Env(5),err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetEmpiricalQuantile(p=0.84_mrk,x=arr,IsXSorted=.true.,q=Env(6),err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetEmpiricalQuantile(p=0.025_mrk,x=arr,IsXSorted=.true.,q=Env(3),err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetEmpiricalQuantile(p=0.975_mrk,x=arr,IsXSorted=.true.,q=Env(4),err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
end subroutine BaRatin_GetEnvelop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_PostProcess(MCMCFile,Nburn,Nread,Nslim,& ! Read properties
                        Hx,& ! list of abscisses Hx where the RC is applied
                        RC,& ! Rating curve object
                        SummaryFile, & ! File containing the statistical summary of MCMC samples
                        SaveSummary,&
                        SpagFile, & ! prefix of File containing all RCs corresponding to MCMC samples
                        EnvFile, & ! prefix of file contaiing envelops extracted from spaghettis
                        SaveSpag,SaveEnvelop,&
                        propagationMatrix,& ! let the user define what is propagated and what is not
                        HQFile,& ! File containing H/Q obs and uncertainties/envelops
                        SaveHQ,&
                        CookedMCMCFile, & ! File containing MCMC samples after burning/slimming/adding b's
                        SaveCooked,&
                        err,mess)! error handling
!^**********************************************************************
!^* Purpose: post-process MCMC samples and write result files for
!^* subsequent use, in particular for subsequent plots
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:16/07/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1. MCMCFile, containing MCMC samples
!^*        2. Nburn, number of lines to discard
!^*        3. Nread, number of lines to read
!^*        4. Nslim, only one line every Nslim will be used
!^*        5. Hx, Hs where the RC will be applied in all post-processings
!^*        6. RC, rating curve object
!^*        7. SummaryFile, Result File containing the statistical summary of MCMC samples
!^*        8. SpagFile, & ! prefix of File containing all RCs corresponding to MCMC samples
!^*        9. EnvFile, & ! prefix of file contaiing envelops extracted from spaghettis
!^*        10. propagationMatrix,& ! let the user define what is propagated and what is not
!^*        11. HQFile, Result File containing H/Q obs and uncertainties/envelops
!^*        12. CookedMCMCFile, File containing MCMC samples after burning/slimming/adding b's
!^* OUT
!^*        1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        2.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile
use utilities_dmsl_kit, only:number_string
use Distribution_tools, only:GAUSS,GetQuantile
use RatingCurve_tools, only:RC_General,ApplyRC

character(*), intent(in)::MCMCFile, SummaryFile, SpagFile,EnvFile,HQFile,CookedMCMCFile
integer(mik), intent(in)::Nburn,Nread,Nslim,propagationMatrix(:,:)
real(mrk), intent(in)::Hx(:)
type(RCType),intent(in):: RC
logical,intent(in)::SaveSpag,SaveEnvelop,SaveSummary,saveHQ,saveCooked
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250), dimension(16),parameter::lefters=(/"N              ",&
                                                   "Minimum        ",&
                                                   "Maximum        ",&
                                                   "Range          ",&
                                                   "Mean           ",&
                                                   "Median         ",&
                                                   "Q10%           ",&
                                                   "Q25%           ",&
                                                   "Q75%           ",&
                                                   "Q90%           ",&
                                                   "St.Dev.        ",&
                                                   "Variance       ",&
                                                   "CV             ",&
                                                   "Skewness       ",&
                                                   "Kurtosis       ",&
                                                   "MaxPost        "/)
character(250), dimension(15),parameter::headers2=(/"Hobs           ",&
                                                    "Hhat(MaxPost)  ",&
                                                    "H2.5%          ",&
                                                    "H97.5%         ",&
                                                    "H16%           ",&
                                                    "H84%           ",&
                                                    "H-sigma        ",&
                                                    "H+sigma        ",&
                                                    "Qobs           ",&
                                                    "Q2.5%          ",&
                                                    "Q97.5%         ",&
                                                    "Q16%           ",&
                                                    "Q84%           ",&
                                                    "Q-sigma        ",&
                                                    "Q+sigma        "/)
character(250),parameter::procname='BaRatin_PostProcess'
integer(mik)::i,j,errcode,k,n,nh,Nmcmc,NHx,p,m,ml(1),nCol
real(mrk), allocatable::Summary(:,:),Env(:,:),DS(:,:)
real(mrk),pointer::mcmc(:,:),SpagTeta(:,:),SpagTot(:,:)
character(250), allocatable::headers(:)
logical::feas
real(mrk)::std,u5,u25,gen
character(250)::fmt

err=0;mess=''

! Read MCMC
call BaRatin_ReadMCMC(trim(MCMCFile),Nburn,Nread,Nslim,& ! Read properties
                        mcmc,err,mess)! error handling
if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
nMCMC=size(mcmc,dim=1);nCol=size(mcmc,dim=2)
ml=maxloc(mcmc(:,RC%nteta+RC%nsigmaf+RC%nHerror+1))
if(saveCooked) then
    open(unit=1,file=trim(CookedMCMCFile),status='replace')
    write(1,*) "Proper headers to be implemented..."
    fmt='('//trim(number_string(nCol))//'e15.6)'
    do i=1,nMCMC;write(1,trim(fmt)) mcmc(i,:);enddo
    close(1)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STEP 1 : SummaryFile for MCMC sampling !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(SaveSummary) then
    if(allocated(Summary)) deallocate(Summary)
    if(allocated(headers)) deallocate(headers)
    if(RC%RCID==RC_General) then
        allocate(Summary(15,RC%nteta+RC%nsigmaf+RC%nHerror+RC%ncontrol-1),headers(RC%nteta+RC%nsigmaf+RC%nHerror+RC%ncontrol-1))
    else
        allocate(Summary(15,RC%nteta+RC%nsigmaf+RC%nHerror),headers(RC%nteta+RC%nsigmaf+RC%nHerror))
    endif
    do i=1,RC%nteta+RC%nsigmaf+RC%nHerror
        call GetEmpiricalStats(x=mcmc(:,i), all15=Summary(:,i),err=err,mess=mess)
        if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    enddo
    if(RC%RCID==RC_General) then
        do i=(RC%nteta+RC%nsigmaf+RC%nHerror+2),(RC%nteta+RC%nsigmaf+RC%nHerror+1+RC%ncontrol-1)
            call GetEmpiricalStats(x=mcmc(:,i), all15=Summary(:,i-1),err=err,mess=mess)
            if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
        enddo
    endif
    ! Create Headers
    do i=1,RC%nteta;headers(i)="Teta"//trim(number_string(i));enddo
    do i=1,RC%nsigmaf;headers(RC%nteta+i)="RemnantSTD"//trim(number_string(i));enddo
    if(RC%Herror) then
        k=0
        do i=1,RC%nobs
            if(RC%Hsigma(i) /= 0._mrk) then
                k=k+1;headers(RC%nteta+RC%nsigmaf+k)="trueH"//trim(number_string(i))
            endif
        enddo
    endif
    if(RC%RCID==RC_General) then
        do i=1,RC%ncontrol-1
            headers(RC%nteta+RC%nsigmaf+RC%nHerror+i)="GeneralRC_b"//trim(number_string(i+1))
        enddo
    endif
    !Write2File
    open(unit=1,file=trim(SummaryFile),status='replace')
    if(RC%RCID==RC_General) then
        nh=RC%nteta+RC%nsigmaf+RC%nHerror+1+RC%ncontrol-1
    else
        nh=RC%nteta+RC%nsigmaf+RC%nHerror+1
    endif
    fmt='('//trim(number_string(nH))//'A15)'
    write(1,trim(fmt)) "Stat           ",headers

    fmt='(A15,'//trim(number_string(nH))//'e15.6)'
    do i=1,15
        write(1,trim(fmt)) lefters(i), Summary(i,:)
    enddo
    if(RC%RCID==RC_General) then
        write(1,trim(fmt)) lefters(16), mcmc(ml(1),1:(RC%nteta+RC%nsigmaf+RC%nHerror)),&
          mcmc(ml(1),(RC%nteta+RC%nsigmaf+RC%nHerror+2):(RC%nteta+RC%nsigmaf+RC%nHerror+1+RC%ncontrol-1))
    else
        write(1,trim(fmt)) lefters(16), mcmc(ml(1),1:(nh-1))
    endif
    close(1)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STEP 2 : spaghetti & envelop files !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,size(propagationMatrix,dim=1)
    call BaRatin_H2Q(RC=RC,& ! Rating curve object
                     Hx=Hx,&
                     mcmc=mcmc,&
                     RemnantOption=0,&
                     PropagateWhat=propagationMatrix(i,:),&
                     SpagFilePrefix=SpagFile,&
                     EnvelopFilePrefix=EnvFile,&
                     SaveSpag=SaveSPag,SaveEnvelop=SaveEnvelop,&
                     err=err,mess=mess)! error handling
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STEP 2 : H-Q data summary !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(SaveHQ) then
    ! get 5% and 25% u-quantiles from a N(0,1)
    ! Replaced by 2.5% and 16% on Jerome request (bloody hydraulicians!)
    call GetQuantile(DistId=GAUSS,p=0.025_mrk,par=(/0._mrk,1._mrk/),q=u5,feas=feas,err=err,mess=mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    call GetQuantile(DistId=GAUSS,p=0.16_mrk,par=(/0._mrk,1._mrk/),q=u25,feas=feas,err=err,mess=mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(DS)) deallocate(DS);allocate(DS(RC%nobs,15))
    k=0
    do i=1,RC%nobs
        DS(i,1)=RC%Hobs(i)
        if(RC%Hsigma(i)==0._mrk) then
            DS(i,2)=RC%Hobs(i)
        else
            k=k+1
            DS(i,2)=mcmc(ml(1),RC%nteta+RC%nsigmaf+k)
        endif
        DS(i,3)=RC%Hobs(i)+RC%Hsigma(i)*u5
        DS(i,4)=RC%Hobs(i)+RC%Hsigma(i)*(-1._mrk*u5)
        DS(i,5)=RC%Hobs(i)+RC%Hsigma(i)*u25
        DS(i,6)=RC%Hobs(i)+RC%Hsigma(i)*(-1._mrk*u25)
        DS(i,7)=RC%Hobs(i)-RC%Hsigma(i)
        DS(i,8)=RC%Hobs(i)+RC%Hsigma(i)
        DS(i,9)=RC%Qobs(i)
        DS(i,10)=RC%Qobs(i)+RC%Qsigma(i)*u5
        DS(i,11)=RC%Qobs(i)+RC%Qsigma(i)*(-1._mrk*u5)
        DS(i,12)=RC%Qobs(i)+RC%Qsigma(i)*u25
        DS(i,13)=RC%Qobs(i)+RC%Qsigma(i)*(-1._mrk*u25)
        DS(i,14)=RC%Qobs(i)-RC%Qsigma(i)
        DS(i,15)=RC%Qobs(i)+RC%Qsigma(i)
    enddo
    open(unit=1,file=trim(HQFile),status='replace')
    write(1,'(15A15)') headers2
    do i=1,RC%nobs;write(1,'(15e15.6)') DS(i,:);enddo
close(1)
endif
end subroutine BaRatin_PostProcess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_Propagation(MCMCFile,Nburn,Nread,Nslim,& ! Read properties
                        Hfile,&                    ! File containing Hs where the RC is applied + Sdev of Hs
                        RC,& ! Rating curve object
                        SpagFile, & ! prefix of File containing all RCs corresponding to MCMC samples
                        EnvFile, & ! prefix of file contaiing envelops extracted from spaghettis
                        SaveSpag,SaveEnvelop,&
                        propagationMatrix,& ! let the user define what is propagated and what is not
                        err,mess)! error handling
!^**********************************************************************
!^* Purpose: Propagates RC uncertainty from H to Q
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:16/07/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1. MCMCFile, containing MCMC samples
!^*        2. Nburn, number of lines to discard
!^*        3. Nread, number of lines to read
!^*        4. Nslim, only one line every Nslim will be used
!^*        5. Hfile, File containing Hs where the RC is applied + Sdev of Hs
!^*        6. RC, rating curve object
!^*        7. SpagFile, & ! prefix of File containing all RCs corresponding to MCMC samples
!^*        8. EnvFile, & ! prefix of file contaiing envelops extracted from spaghettis
!^*        9. propagationMatrix,& ! let the user define what is propagated and what is not
!^* OUT
!^*        1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        2.mess, error message
!^**********************************************************************
use DataRW_tools, only:DatRead
use utilities_dmsl_kit,only:number_string
character(*), intent(in)::MCMCFile, HFile, SpagFile,EnvFile
integer(mik), intent(in)::Nburn,Nread,Nslim,propagationMatrix(:,:)
type(RCType),intent(in):: RC
logical,intent(in)::SaveSpag,SaveEnvelop
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaRatin_Propagation'
integer(mik),parameter::ncol=4
integer(mik)::Nmcmc,NHx,i
real(mrk),pointer::mcmc(:,:),HfileContent(:,:)
character(250),dimension(ncol)::headers

err=0;mess=''
! Read H-file
call DatRead(file=Hfile,ncol=ncol,y=HfileContent,headers=headers,err=err,mess=mess)
if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif

! Read MCMC
call BaRatin_ReadMCMC(trim(MCMCFile),Nburn,Nread,Nslim,& ! Read properties
                        mcmc,err,mess)! error handling
if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
nMCMC=size(mcmc,dim=1)

! spaghettis & envelops
do i=1,size(propagationMatrix,dim=1)
    write(*,'(A50)') 'Propagation experiment: '//trim(number_string(i))//&
                     &' / '//trim(number_string(size(propagationMatrix,dim=1)))
    call BaRatin_H2Q(RC=RC,& ! Rating curve object
                        Hx=HfileContent(:,1),&
                        Hx_Std=HfileContent(:,2),&
                        Hx_Bias=HfileContent(:,4),&
                        Hx_BiasIndx=nint(HfileContent(:,3)),&
                        mcmc=mcmc,&
                        PropagateWhat=propagationMatrix(i,:),&
                        PrintCounter=.true.,&
                        SpagFilePrefix=SpagFile,&
                        EnvelopFilePrefix=EnvFile,&
                        SaveSpag=SaveSpag,SaveEnvelop=SaveEnvelop,&
                        err=err,mess=mess)! error handling
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    write(*,'(A50)') 'ooooooooooooooooooooooooooooooooooooooooooo'
enddo
end subroutine BaRatin_Propagation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine Sigmafunk_GetParNumber(funk, npar, err, mess)

!^**********************************************************************
!^* Purpose: Get number of parameters of the selected Sigmafunk
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 29/04/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.funk, which function?? (e.g., 'Linear')
!^* OUT
!^*        1.npar
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
character(*), intent(in)::funk
integer(mik), intent(out)::npar
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Sigmafunk_GetParNumber'

err=0;mess='';npar=undefIN
select case(trim(funk))
case('Constant','Proportional')
    npar=1
case('Linear')
    npar=2
case('Exponential', 'Gaussian')
    npar=3
case default
    err=1;mess=trim(procname)//':'//trim(Baratin_message(6))
end select

end subroutine Sigmafunk_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine Sigmafunk_Apply(funk, par, Qrc, res, err, mess)

!^**********************************************************************
!^* Purpose: Apply the selected Sigmafunk
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 29/04/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.funk, which function?? (e.g., 'Constant','Linear')
!^*        2.par, parameters of funk
!^*        3.Qrc, Q given by the rating curve
!^* OUT
!^*        1.res, result
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
character(*), intent(in)::funk
real(mrk), intent(in)::par(:),Qrc
real(mrk), intent(out)::res
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Sigmafunk_Apply'
integer(mik)::npar

err=0;mess='';res=undefRN
call Sigmafunk_GetParNumber(funk, npar, err, mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(par)/=npar) then;err=1;mess=trim(procname)//trim(Baratin_message(5));return;endif

select case(trim(funk))
case('Constant')
    res=par(1)
case('Proportional')
    res=par(1)*Qrc
case('Linear')
    res=par(1)+par(2)*Qrc
case('Exponential')
    res=par(1)+( par(3)-par(1) )*( 1._mrk - exp(- (Qrc/par(2))**1 ) )
case('Gaussian')
    res=par(1)+( par(3)-par(1) )*( 1._mrk - exp(- (Qrc/par(2))**2 ) )
case default
    err=1;mess=trim(procname)//trim(Baratin_message(6))
end select

end subroutine SigmaFunk_Apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_Workspace(file,&
                                   workspace,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_BaRatin
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_RunOptions
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::workspace,mess
! locals
character(250),parameter::procname='BaRatinConfig_Read_Workspace'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) workspace
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine BaRatinConfig_Read_Workspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_Data(file,&
                                   DataFile,nHeader,nobs,ncol,HobsCol,HsigmaCol,&
                                   QobsCol,QsigmaCol,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_Data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 26/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_Data
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err,nHeader,nobs,ncol,HobsCol,HsigmaCol,QobsCol,QsigmaCol
character(*),intent(out)::mess,DataFile
! locals
character(250),parameter::procname='BaRatinConfig_Read_Data'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) DataFile
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nHeader
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nobs
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) ncol
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) HobsCol
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) HsigmaCol
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) QobsCol
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) QsigmaCol
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine BaRatinConfig_Read_Data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_MCMC(file,teta0,RemnantSigma0,smallStd,&
                                   nAdapt,nCycles,Nburn,Nread,nSlim,InitStdMode,&
                                   BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                                   teta_std0,RemnantSigma_std0,&
                                   saveMCMC,MCMCfile,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file,teta0,RemnantSigma0,smallStd
!^* OUT
!^*        1. All parameters in Config_MCMC + starting stds
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
real(mrk), intent(in)::teta0(:),RemnantSigma0(:),smallStd
integer(mik), intent(out)::err,nAdapt,nCycles,Nburn,Nread,nSlim,InitStdMode
real(mrk),intent(out)::BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult
real(mrk), pointer:: teta_std0(:),RemnantSigma_std0(:) ! starting point/std for teta
character(*),intent(out)::mess,MCMCfile
logical, intent(out)::saveMCMC
! locals
character(250),parameter::procname='BaRatinConfig_Read_MCMC'
integer(mik)::unt,nteta,i
real(mrk),allocatable::stdFactorV(:)
real(mrk)::stdFactor

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
! general properties of the sampler
read(unt,*,iostat=err) nAdapt
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nCycles
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) BurnFactor
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
Nburn=(nAdapt*nCycles)*BurnFactor
Nread=nAdapt*nCycles-Nburn
read(unt,*,iostat=err) nSlim
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) MinMoveRate
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) MaxMoveRate
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DownMult
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) UpMult
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) InitStdMode
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
! Handling of initial jump stdev
nteta=size(teta0)
if(associated(teta_std0)) nullify(teta_std0);allocate(teta_std0(nTeta))
if(associated(RemnantSigma_std0)) nullify(RemnantSigma_std0);allocate(RemnantSigma_std0(size(RemnantSigma0)))
if(InitStdMode==1) then
    if(allocated(stdFactorV)) deallocate(stdFactorV);allocate(stdFactorV(nTeta))
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) stdFactorV
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    FORALL (i=1:nTeta) teta_std0(i) = stdFactorV(i)*abs(teta0(i))
    where(teta_std0==0._mrk) teta_std0=smallStd
    if(allocated(stdFactorV)) deallocate(stdFactorV);allocate(stdFactorV(size(RemnantSigma0)))
    read(unt,*,iostat=err) stdFactorV
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    FORALL (i=1:size(RemnantSigma0)) RemnantSigma_std0(i) = stdFactorV(i)*abs(RemnantSigma0(i))
    where(RemnantSigma_std0==0._mrk) RemnantSigma_std0=smallStd
else
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) stdFactor
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    teta_std0=stdFactor*abs(teta0)
    where(teta_std0==0._mrk) teta_std0=smallStd
    RemnantSigma_std0 = stdFactor*abs(RemnantSigma0)
    where(RemnantSigma_std0==0._mrk) RemnantSigma_std0=smallStd
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
endif
! result file
read(unt,*,iostat=err) !cosmetics
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) saveMCMC
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) MCMCfile
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif

close(unt)

end subroutine BaRatinConfig_Read_MCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_RatingCurve(file,&
                                   RCID,nTeta,&
                                   teta0,&
                                   PriorList_teta,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_RatingCurve
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_RatingCurve + prior lists
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use Distribution_tools, only:GetParNumber
character(*), intent(in)::file
integer(mik), intent(out)::err,nTeta
real(mrk), pointer::teta0(:)
Type(PriorListType), pointer::PriorList_teta(:)
character(*),intent(out)::mess,RCID
! locals
character(250),parameter::procname='BaRatinConfig_Read_RatingCurve'
integer(mik)::unt,i,np
logical::exist

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) RCID
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nTeta
if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
if(associated(PriorList_teta)) nullify(PriorList_teta);allocate(PriorList_teta(nTeta))
if(associated(teta0)) nullify(teta0);allocate(teta0(nTeta))
do i=1,nteta
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) teta0(i)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) PriorList_teta(i)%dist
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    PriorList_teta(i)%dist=trim(PriorList_teta(i)%dist)
    call GetParNumber(DistID=PriorList_teta(i)%dist, npar=np, err=err, mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(PriorList_teta(i)%par)) deallocate(PriorList_teta(i)%par)
    allocate(PriorList_teta(i)%par(np))
    read(unt,*,iostat=err) PriorList_teta(i)%par(:)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)

end subroutine BaRatinConfig_Read_RatingCurve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_ControlMatrix(file,nTeta,nElementalRC,nHeaderControlMatrix,&
                                   ControlMatrix,ncontrol,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read BaRatinConfig_Read_ControlMatrix
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_ControlMatrix
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit,number_string
use RatingCurve_tools, only:RC_General_CheckControlMatrix
character(*), intent(in)::file
integer(mik), intent(in)::nTeta,nElementalRC,nHeaderControlMatrix
integer(mik),pointer::ControlMatrix(:,:)
integer(mik), intent(out)::err,ncontrol
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatinConfig_Read_ControlMatrix'
integer(mik)::unt
logical::feas

err=0;mess=''

if( MOD(nteta,nElementalRC)/=0 ) then
    err=1
    mess=trim(procname)//':'//trim(BaRatin_Message(1))//trim(number_string(nElementalRC))
endif

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif

ncontrol=nTeta/nElementalRC
if(associated(ControlMatrix)) nullify(ControlMatrix);allocate(ControlMatrix(ncontrol,ncontrol))
call ReadData_i(DataFile=file,unt=unt,&
                nrow=ncontrol,ncol=ncontrol,nHeader=nHeaderControlMatrix,&
                X=ControlMatrix,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
call RC_General_CheckControlMatrix(ControlMatrix,feas,err,mess)
if((.not.feas).or.(err/=0)) then;mess=trim(procname)//':'//trim(mess);return;endif
close(unt)
end subroutine BaRatinConfig_Read_ControlMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_RemnantSigma(file,RCID,ncontrol,&
                                   RemnantSigma_funk,&
                                   RemnantSigma0,PriorList_RemnantSigma,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_RemnantSigma
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file,RCID,ncontrol,
!^* OUT
!^*        1. All parameters in Config_ControlMatrix + prior lists
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use Distribution_tools, only:GetParNumber
use RatingCurve_tools, only:RC_General
character(*), intent(in)::file,RCID
integer(mik), intent(in)::ncontrol
Type(PriorListType), pointer:: PriorList_RemnantSigma(:)
real(mrk),pointer:: RemnantSigma0(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess,RemnantSigma_funk
! locals
character(250),parameter::procname='BaRatinConfig_Read_RemnantSigma'
integer(mik)::unt,i,SigmaFunk_np,np,SigmaFunk_np0
logical::exist

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) RemnantSigma_funk
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SigmaFunk_np
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
! check number of parameters is correct
call Sigmafunk_GetParNumber(funk=RemnantSigma_funk, npar=SigmaFunk_np0, err=err, mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(SigmaFunk_np0/=SigmaFunk_np) then
    err=1;mess=trim(procname)//':'//Baratin_message(4);return
endif
if(associated(RemnantSigma0)) nullify(RemnantSigma0);allocate(RemnantSigma0(SigmaFunk_np))
if(associated(PriorList_RemnantSigma)) nullify(PriorList_RemnantSigma);allocate(PriorList_RemnantSigma(SigmaFunk_np))
!read each parameter block
do i=1,SigmaFunk_np
    read(unt,*,iostat=err)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) RemnantSigma0(i)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) PriorList_RemnantSigma(i)%dist
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    call GetParNumber(DistID=PriorList_RemnantSigma(i)%dist, npar=np, err=err, mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(PriorList_RemnantSigma(i)%par)) deallocate(PriorList_RemnantSigma(i)%par)
    allocate(PriorList_RemnantSigma(i)%par(np))
    read(unt,*,iostat=err) PriorList_RemnantSigma(i)%par(:)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)

end subroutine BaRatinConfig_Read_RemnantSigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_RunOptions(file,&
                                   DoMCMC,DoPostProcess,DoPriorRC,DoH2Q,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_RunOptions
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_RunOptions
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
logical, intent(out)::DoMCMC,DoPostProcess,DoPriorRC,DoH2Q
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatinConfig_Read_RunOptions'
integer(mik)::unt
logical::exist

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
inquire(file=trim(file),exist=exist)
if(exist) then !
    open(unit=unt,file=trim(file), status='old', iostat=err)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
    read(unt,*,iostat=err) DoPriorRC
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoMCMC
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoPostProcess
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoH2Q
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    close(unt)
else ! default settings - MCMC+PostProcess only
    DoMCMC=.true.
    DoPostProcess=.true.
    DoPriorRC=.false.
    DoH2Q=.false.
endif

end subroutine BaRatinConfig_Read_RunOptions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_PriorRC(file,Nsim,Hx,&
                                      SaveSpag,SpagPrefix,&
                                      SaveEnvelop,EnvelopPrefix,&
                                      err,mess)

!^**********************************************************************
!^* Purpose: Read Config_PriorRC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 17/07/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_PriorRC
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::Nsim,err
character(*),intent(out)::SpagPrefix,EnvelopPrefix,mess
real(mrk),pointer::Hx(:)
logical,intent(out)::SaveSpag,SaveEnvelop
! locals
character(250),parameter::procname='BaRatinConfig_Read_PriorRC'
integer(mik)::unt,nHx,i
real(mrk)::Hmin,Hmax,Hstep

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) Nsim
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Hmin,Hmax,Hstep
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
nHx=ceiling((Hmax-Hmin)/Hstep)+1
if(associated(Hx)) nullify(Hx);allocate(Hx(nHx))
Hx=(/(Hmin+i*Hstep,i=0,NHx-1)/)
read(unt,*,iostat=err) saveSpag
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SpagPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) saveEnvelop
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) EnvelopPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine BaRatinConfig_Read_PriorRC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_PostProcessing(file,&
                                   SaveCooked,CookedFile,&
                                   SaveSummary,SummaryFile,&
                                   SaveHQ,HQfile,&
                                   Hx,&
                                   SaveSpag,SpagPrefix,&
                                   SaveEnvelop,EnvelopPrefix,&
                                   PropMatrix,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_PostProcessing
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 17/07/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_PostProcessing
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::SpagPrefix,EnvelopPrefix,SummaryFile,CookedFile,HQfile,mess
real(mrk),pointer::Hx(:)
integer(mik),pointer::propMatrix(:,:)
logical,intent(out)::SaveCooked,SaveSpag,SaveEnvelop,SaveSummary,SaveHQ
! locals
character(250),parameter::procname='BaRatinConfig_Read_PostProcessing'
integer(mik),parameter::ncol=3
integer(mik)::unt,nHx,i,nrow
real(mrk)::Hmin,Hmax,Hstep

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
! MCMC cooking
read(unt,*,iostat=err) !cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SaveCooked
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) CookedFile
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
! MCMC summary
read(unt,*,iostat=err) !cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SaveSummary
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SummaryFile
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
! HQ summary
read(unt,*,iostat=err) !cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SaveHQ
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) HQFile
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
! Posterior rating curves
read(unt,*,iostat=err) !cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Hmin,Hmax,Hstep
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
nHx=ceiling((Hmax-Hmin)/Hstep)+1
if(associated(Hx)) nullify(Hx);allocate(Hx(nHx))
Hx=(/(Hmin+i*Hstep,i=0,NHx-1)/)
read(unt,*,iostat=err) saveSpag
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SpagPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) saveEnvelop
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) EnvelopPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nrow
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
if(associated(PropMatrix)) nullify(PropMatrix);allocate(PropMatrix(nrow,ncol))
PropMatrix=0
do i=1,nrow
    read(unt,*,iostat=err) PropMatrix(i,(/tcol,rcol/))
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)

end subroutine BaRatinConfig_Read_PostProcessing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_H2QPropagation(file,&
                                   HHfile,&
                                   SaveSpag,SpagPrefix,&
                                   SaveEnvelop,EnvelopPrefix,&
                                   PropMatrix,&
                                   err,mess)

!^**********************************************************************
!^* Purpose: Read Config_H2QPropagation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 17/07/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_H2QPropagation
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::HHfile,SpagPrefix,EnvelopPrefix,mess
integer(mik),pointer::propMatrix(:,:)
logical,intent(out)::SaveSpag,SaveEnvelop
! locals
character(250),parameter::procname='BaRatinConfig_Read_H2QPropagation'
integer(mik),parameter::ncol=3
integer(mik)::unt,i,nrow

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif

read(unt,*,iostat=err) HHfile
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) saveSpag
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) SpagPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) saveEnvelop
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) EnvelopPrefix
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nrow
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
if(associated(PropMatrix)) nullify(PropMatrix);allocate(PropMatrix(nrow,ncol))
do i=1,nrow
    read(unt,*,iostat=err) PropMatrix(i,:)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)

end subroutine BaRatinConfig_Read_H2QPropagation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatinConfig_Read_Aggregation(file,&
                                   Tfile_IN,Tfile_nhead,nT,Tformat_IN,Tfile_OUT,Tformat_OUT,&
                                   Sfile_IN,Sfile_nhead,nspag,mvcode,Sfile_OUT,&
                                   FirstDate,LastDate,scale,scaleMult,&
                                   funk,safeMV,safeInterpol,option,&
                                   err,mess)
!^**********************************************************************
!^* Purpose: Read Config_BaRatinAggregator
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 14/05/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file
!^* OUT
!^*        1. All parameters in Config_BaRatinAggregator
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
real(mrk),intent(out)::mvcode
integer(mik), intent(out)::Tfile_nhead,Sfile_nhead,nT,Tformat_IN,Tformat_OUT,nspag,&
                           firstDate(6),LastDate(6),scaleMult,option,err
character(*),intent(out)::Tfile_IN,Tfile_OUT,Sfile_IN,Sfile_OUT,scale,&
                          funk,mess
logical,intent(out)::safeMV,safeInterpol
! locals
character(250),parameter::procname='BaRatinConfig_Read_Aggregation'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif

read(unt,*,iostat=err) ! Cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Tfile_IN
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Tfile_nhead
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nT
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Tformat_IN
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Tfile_OUT
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Tformat_OUT
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Sfile_IN
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Sfile_nhead
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nspag
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) mvcode
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Sfile_OUT
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) ! Cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) FirstDate
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) LastDate
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) scale
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) scalemult
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) ! Cosmetics
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) funk
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) safeMV
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) safeInterpol
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) option
if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
close(unt)
end subroutine BaRatinConfig_Read_Aggregation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function BaRatin_Message(id)
! returns error/warning messages
integer(mik),intent(in)::id
character(250)::BaRatin_Message

select case(id)
case(1)
    BaRatin_Message='number of parameters for a RC of type RC_General has to be a multiple of:'
case(2)
    BaRatin_Message='problem reading MCMC file during burnin'
case(3)
    BaRatin_Message='WARNING:problem reading MCMC file before [Nread] lines'
case(4)
    BaRatin_Message='Incorrect number of parameters for remnant sigma function'
case(5)
    BaRatin_Message='size mismacth [funk,par]'
case(6)
    BaRatin_Message='unknown function [funk]'
case(7)
    BaRatin_Message='remnant-propagation is not allowed in prior mode'
case(8)
    BaRatin_Message='[Hx_Std] compulsory when h-uncertainty is propagated'
case(9)
    BaRatin_Message='[nsim] compulsory when [mcmc] is not provided'
case(10)
    BaRatin_Message='prior mode is not feasible'
case(11)
    BaRatin_Message='impossible to generate from prior'
case(12)
    BaRatin_Message='size mismatch [teta]'
case(13)
    BaRatin_Message='size mismatch [RemnantSigma]'
case(14)
    BaRatin_Message='Fatal: size mismatch [x]'
case(15)
    BaRatin_Message='size mismatch in input stage data'
case(16)
    BaRatin_Message='size mismatch [PriorList_RemnantSigma]'
case(17)
    BaRatin_Message='size mismatch [PriorList_teta]'
case(18)
    BaRatin_Message='size mismatch [teta / RCID]'
case(19)
    BaRatin_Message='size mismatch [RemnantSigma / RemnantSigma_funk]'
case(20)
    BaRatin_Message='size mismatch in input stage or streamflow data'
case default
    BaRatin_Message='unknown message ID'
end select

end function BaRatin_Message

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_ConsoleMessage(id,mess)
! print the the desired error message on screen
integer(mik),intent(in)::id
character(*), intent(in)::mess

select case(id)
case (1) ! Welcome!
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) '**           BaRatin           **'
    write(*,*) '*********************************'
    write(*,*) '* Bayesian rating curve fitting *'
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) ''
    write(*,*) 'Please enter the workspace folder'
    write(*,*) 'containing configuration files:'
case(2)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**    Prior Rating Curve...    **'
    write(*,*) '*********************************'
    write(*,*) ''
case(3)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**  Prior Rating Curve : DONE! **'
    write(*,*) '*********************************'
    write(*,*) ''
case(4)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**     MCMC sampling...        **'
    write(*,*) '*********************************'
    write(*,*) ''
case (5)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**    MCMC sampling : DONE!    **'
    write(*,*) '*********************************'
    write(*,*) ''
    write(*,*) 'Open MCMC file and monitor results...'
    write(*,*) 'MCMC File = ', trim(mess)
case(6)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**     Post-processing...      **'
    write(*,*) '*********************************'
    write(*,*) ''
case(7)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**   Post-processing : DONE!   **'
    write(*,*) '*********************************'
    write(*,*) ''
case(8)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**     H2Q propagation...      **'
    write(*,*) '*********************************'
    write(*,*) ''
case(9)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**   H2Q propagation : DONE!   **'
    write(*,*) '*********************************'
    write(*,*) ''
case(10)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**    Reading config files..   **'
    write(*,*) '*********************************'
    write(*,*) ''
case(11)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '** Reading config files: DONE! **'
    write(*,*) '*********************************'
    write(*,*) ''
case (101) ! Aggregator welcome
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) '**    BaRatin - Aggregator     **'
    write(*,*) '*********************************'
    write(*,*) '*   Aggregate your spaghettis!  *'
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) ''
case(999) ! Ciao!
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) '**           BaRatin           **'
    write(*,*) '*********************************'
    write(*,*) '*           All done!           *'
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) ''
    write(*,*) 'Thanx for using me...'
    write(*,*) 'Press [enter] and I''ll go away'
case (-1) ! Fatal Error - general
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-2) ! Fatal Error - Open File
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while opening the following file."
    write(*,*) trim(mess)
    write(*,*) "Please check path to file."
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-3) ! Fatal Error - reading a Config file
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while reading the config file:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-4) ! Fatal Error - Prior RC
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while generating the prior rating curve."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-5) ! Fatal Error - MCMC fit
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while fitting the rating curve."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-6) ! Fatal Error - Post-processing
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while post-processing MCMC samples."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-7) ! Fatal Error - propagation
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while propagating uncertainty to hydrograph"
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case(-8) ! Fatal Error - writing to a file
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "while writting to the file:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
case default
    write(*,*) ""
    write(*,*) "BaRatin: a FATAL ERROR has occured"
    write(*,*) "with unknown error ID..."
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaRatin_Fatal_Exit
end select

end subroutine BaRatin_ConsoleMessage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_Fatal_Exit
!read(*,*)
!STOP
call exit(1)
end subroutine BaRatin_Fatal_Exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaRatin_ReadData(file,nrow,ncol,nHeader,&
                            HobsCol,HsigmaCol,QobsCol,QsigmaCol,&
                            Hobs,Qobs,Hsigma,Qsigma,err,mess)

!^**********************************************************************
!^* Purpose: Read data file
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 26/09/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.file, DataFile
!^*        2.nrow (WITHOUT header lines)
!^*        3.ncol
!^*        4.nHeader, number of skipped header lines
!^* OUT
!^*        1.Hobs,Qobs,Hsigma,Qsigma
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(in)::nrow,ncol,nHeader,HobsCol,HsigmaCol,QobsCol,QsigmaCol
real(mrk),pointer::Hobs(:),Qobs(:),Hsigma(:),Qsigma(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaRatin_ReadData'
integer(mik)::i,unt
real(mrk)::X(nrow,ncol)

err=0;mess='';

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

open(unit=unt,file=trim(file), status='old',iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif

do i=1,nHeader
    read(unt,*,iostat=err)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo

do i=1,nrow
    read(unt,*,iostat=err) X(i,:)
    if(err/=0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)
if(associated(Hobs)) nullify(Hobs);
if(associated(Qobs)) nullify(Qobs);
if(associated(Hsigma)) nullify(Hsigma);
if(associated(Qsigma)) nullify(Qsigma);
allocate(Hobs(nrow),Qobs(nrow),Hsigma(nrow),Qsigma(nrow))
Hobs=X(:,HobsCol)
Qobs=X(:,QobsCol)
Hsigma=X(:,HsigmaCol)
Qsigma=X(:,QsigmaCol)
end subroutine BaRatin_ReadData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadData_i(DataFile,unt,nrow,ncol,nHeader,X,err,mess)

!^**********************************************************************
!^* Purpose: Read matrix-like data file
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:04/10/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1.DataFile
!^*        2.nrow (WITHOUT header lines)
!^*        3.ncol
!^*        4.nHeader, number of skipped header lines
!^* OUT
!^*        1.X
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
character(*), intent(in)::DataFile
integer(mik), intent(in)::nrow,ncol,nHeader,unt
integer(mik), intent(out)::X(nrow,ncol)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='ReadData_i'
integer(mik)::i

err=0;mess='';X=undefIN

do i=1,nHeader
    read(unt,*,iostat=err)
    if(err/=0) then;call Baratin_ConsoleMessage(messID_Read,trim(DataFile));endif
enddo

do i=1,nrow
    read(unt,*,iostat=err) X(i,:)
    if(err/=0) then;call Baratin_ConsoleMessage(messID_Read,trim(DataFile));endif
enddo

end subroutine ReadData_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module BaRatin_tools

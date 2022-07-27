program main

! use what?
use EmpiricalStats_tools
use kinds_dmsl_kit
use RatingCurve_tools
use BaRatin_tools
use BayesianEstimation_tools, only:PriorListType
use Distribution_tools
use numerix_dmsl_kit, only: quicksort
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constants
real(mrk), parameter:: smallStd=0.1_mrk
integer(mik), parameter:: nElementalRC=3,& ! number of parameters of the elemental RC equation (power)
                          nHeaderControlMatrix=0 ! header lines for the Control Matrix File
character(250),parameter::  Config_BaRatin="Config_BaRatin.txt",&
                            Config_data="Config_Data.txt",&
                            Config_RC="Config_RatingCurve.txt",&
                            Config_MCMC="Config_MCMC.txt",&
                            Config_ControlMatrix="Config_ControlMatrix.txt",&
                            Config_RemnantSigma="Config_RemnantSigma.txt",&
                            Config_RunOptions="Config_RunOptions.txt",&
                            Config_PriorRC="Config_PriorRC.txt",&
                            Config_PostProcessing="Config_PostProcessing.txt",&
                            Config_H2QPropagation="Config_H2QPropagation.txt"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BaRatin tuning - Variables to be read from input file
character(250):: DataFile ! path to datafile
character(250):: HHFile ! path to limnograph
integer(mik):: nobs       ! number of obs in Data file
integer(mik):: ncol       ! number of columns in Data file
integer(mik):: nHeader    ! number of header lines to skip in the data file
character(250):: MCMCfile ! path to file where MCMC samples will be stored
character(250):: SummaryFile !path for MCMC summary
character(250):: CookedFile !path for cooked MCMC samples
character(250):: HQFile ! path for HQ summary
integer(mik):: HobsCol,HsigmaCol,QobsCol,QsigmaCol ! columns for Hobs,Hsigma,Qobs,Qsigma in DataFile
character(100):: RCID !ID of the rating curve
integer(mik)::nTeta ! number of parameters for the chosen RC
Type(PriorListType), pointer::PriorList_teta(:)
Type(PriorListType), pointer:: PriorList_RemnantSigma(:)
logical::DoMCMC,DoPostProcess,DoPriorRC,DoH2Q,saveMCMC,saveSpag_prior,&
         saveEnvelop_prior,SaveSummary,SaveHQ,saveSpag_post,saveEnvelop_post,&
         saveSpag_prop,saveEnvelop_prop,saveCooked
integer(mik)::PriorNsim ! Nsim for generating prior RC
real(mrk),pointer::Hx(:),Hx_prior(:)
integer(mik),pointer::propmatrix_post(:,:),propmatrix_prop(:,:)
character(250)::SpagPrefix_prior,EnvelopPrefix_prior,SpagPrefix_post,EnvelopPrefix_post,&
                SpagPrefix_prop,EnvelopPrefix_prop
! MCMC tuning
real(mrk), pointer:: teta0(:),teta_std0(:) ! starting point/std for teta
real(mrk),pointer:: RemnantSigma0(:), RemnantSigma_std0(:)   ! starting point/std for remnant std
integer(mik):: nAdapt,nCycles,Nburn,Nread,Nslim,InitStdMode
real(mrk):: MinMoveRate,MaxMoveRate,DownMult,UpMult, BurnFactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! locals
integer(mik)::i, err, nHerror,ncontrol
integer(mik),pointer::ControlMatrix(:,:)
character(500)::mess, workspace,RemnantSigma_funk
character(1)::dirsep
real(mrk),pointer::Hobs(:),Hsigma(:),Qobs(:),Qsigma(:)
logical(mlk)::feas,exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! determine dirsep according to OS
!DEC$ IF DEFINED(_WIN32)
    dirsep=dirSepCH_win
!DEC$ ELSEIF DEFINED(__linux)
    dirsep=dirSepCH_lix
!DEC$ ELSE
    dirsep=dirSepCH_win
!DEC$ ENDIF

! Welcome user and read Workspace
call BaRatin_ConsoleMessage(1,'')
call BaRatinConfig_Read_Workspace(trim(Config_BaRatin),workspace,err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
!read(*,'(A250)') workspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read compulsory config files
call BaRatin_ConsoleMessage(10,'')

! RunOptions configuration file (optional, see defaults in sub)
call BaRatinConfig_Read_RunOptions(trim(workspace)//dirsep//trim(Config_RunOptions),&
                                   DoMCMC,DoPostProcess,DoPriorRC,DoH2Q,&
                                   err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! DATA configuration file
call BaRatinConfig_Read_Data(trim(workspace)//dirsep//trim(Config_data),&
                             DataFile,nHeader,nobs,ncol,HobsCol,HsigmaCol,&
                             QobsCol,QsigmaCol,&
                             err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! Read Data
call BaRatin_ReadData(DataFile,nobs,ncol,nHeader,&
                      HobsCol,HsigmaCol,QobsCol,QsigmaCol,&
                      Hobs,Qobs,Hsigma,Qsigma,err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! RATING CURVE configuration file
call BaRatinConfig_Read_RatingCurve(trim(workspace)//dirsep//trim(Config_RC),&
                                   RCID,nTeta,&
                                   teta0,&
                                   PriorList_teta,&
                                   err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! RC_General ONLY: read Control Matrix from file
if(RCID==RC_General) then
    call BaRatinConfig_Read_ControlMatrix(trim(workspace)//dirsep//trim(Config_ControlMatrix),&
                                          nTeta,nElementalRC,nHeaderControlMatrix,&
                                          ControlMatrix,ncontrol,&
                                          err,mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
endif

! Config_RemnantSigma
call BaRatinConfig_Read_RemnantSigma(file=trim(workspace)//dirsep//trim(Config_RemnantSigma),&
                                   RCID=RCID,ncontrol=ncontrol,&
                                   RemnantSigma_funk=RemnantSigma_funk,&
                                   RemnantSigma0=RemnantSigma0,&
                                   PriorList_RemnantSigma=PriorList_RemnantSigma,&
                                   err=err,mess=mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! MCMC Configuration file
call BaRatinConfig_Read_MCMC(trim(workspace)//dirsep//trim(Config_MCMC),teta0,RemnantSigma0,smallStd,&
                             nAdapt,nCycles,Nburn,Nread,nSlim,InitStdMode,&
                             BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                             teta_std0,RemnantSigma_std0,&
                             saveMCMC,MCMCfile,&
                             err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

call BaRatin_ConsoleMessage(11,'')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! load rating curve object
call BaRatin_LoadRCobject(Hobs,Hsigma,Qobs,Qsigma,& ! observations and uncertainties
                        RCID,&               ! Rating Curve ID
                        ControlMatrix,& ! optional Control Matrix for Laurent's general formalisation
                        RemnantSigma_funk,& ! [option 2] chosen f in sdev=f(Qrc)
				        PriorList_teta,PriorList_RemnantSigma,&
				        err,mess)! error handling
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(DoPriorRC) then
    call BaRatin_ConsoleMessage(2,'')
    call BaRatinConfig_Read_PriorRC(file=trim(workspace)//dirsep//trim(Config_PriorRC),&
                                   Nsim=PriorNsim,Hx=Hx_prior,&
                                   SaveSpag=SaveSpag_prior,SpagPrefix=SpagPrefix_prior,&
                                   SaveEnvelop=SaveEnvelop_prior,EnvelopPrefix=EnvelopPrefix_prior,&
                                   err=err,mess=mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
    call BaRatin_H2Q(RC=RC,& ! Rating curve object
                        Hx=Hx_prior,&
                        nsim=PriorNsim,&
                        PropagateWhat=(/0,1,0/),&
                        SpagFilePrefix=trim(workspace)//dirsep//trim(SpagPrefix_prior),&
                        EnvelopFilePrefix=trim(workspace)//dirsep//trim(EnvelopPrefix_prior),&
				        SaveSpag=SaveSpag_prior,&
				        SaveEnvelop=SaveEnvelop_prior,&
				        err=err,mess=mess)! error handling
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
    call BaRatin_ConsoleMessage(3,'')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (DoMCMC) then
    ! MCMC sampling
    call BaRatin_ConsoleMessage(4,'')
    nHerror=count(Hsigma/=0._mrk)
    write(*,*) "Data File = ", trim(DataFile)
    write(*,*) "nobs = ", nobs
    write(*,*) "Rating curve = ", trim(RCID)
    write(*,*) "Number of parameters = ", nteta
    if (nHerror == 0) then
        write(*,*) "All stage values assumed error-free"
    else
        write(*,*) "Number of non-error-free stage values: ", nHerror
    endif
    write(*,*) "Total number of MCMC iterations", nAdapt*nCycles
    write(*,*) "MCMC is running, be patient..."
    call BaRatin_Fit(RC=RC,&
                     !!!!!!! Tuning of the MCMC sampler !!!!!!!!
                     teta0=teta0,RemnantSigma0=RemnantSigma0, &!initial values for teta and remnant std
				     teta_std0=teta_std0,RemnantSigma_std0=RemnantSigma_std0,& ! initial values for the std of jump distribution
				     nAdapt=nAdapt,nCycles=nCycles,&
				     MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				     DownMult=DownMult,UpMult=UpMult,&
                     !!!!!!! END Tuning of the MCMC sampler !!!!!!!!
				     OutFile=trim(workspace)//dirsep//trim(MCMCfile), & ! Output file (for MCMC samples)
				     err=err,mess=mess)! error handling
    if(err>0) then; call BaRatin_ConsoleMessage(-5,trim(mess));endif
    if(.not.saveMCMC) then
        open(unit=666,file=trim(workspace)//dirsep//trim(MCMCfile),status='old')
        close(unit=666,status='delete')
    endif
    call BaRatin_ConsoleMessage(5, trim(trim(workspace)//dirsep//trim(MCMCfile)))
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(DoPostProcess) then
    ! Post-Processing
    call BaRatin_ConsoleMessage(6,'')
    call BaRatinConfig_Read_PostProcessing(file=trim(workspace)//dirsep//trim(Config_PostProcessing),&
                                   SaveCooked=SaveCooked,CookedFile=CookedFile,&
                                   SaveSummary=SaveSummary,SummaryFile=SummaryFile,&
                                   SaveHQ=SaveHQ,HQfile=HQfile,&
                                   Hx=Hx,&
                                   SaveSpag=SaveSpag_post,SpagPrefix=SpagPrefix_post,&
                                   SaveEnvelop=SaveEnvelop_post,EnvelopPrefix=EnvelopPrefix_post,&
                                   PropMatrix=PropMatrix_post,&
                                   err=err,mess=mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-3,trim(mess));endif
    call BaRatin_PostProcess(MCMCFile=trim(workspace)//dirsep//trim(MCMCfile),&
                            Nburn=Nburn,Nread=Nread,Nslim=Nslim,& ! Read properties
                            Hx=Hx,&
                            RC=RC,&
                            SummaryFile=trim(workspace)//dirsep//trim(SummaryFile), & ! File containing the statistical summary of MCMC samples
                            SaveSummary=SaveSummary,&
                            SpagFile=trim(workspace)//dirsep//trim(SpagPrefix_post), & ! File containing all RCs corresponding to MCMC samples - teta uncertainty only
                            EnvFile=trim(workspace)//dirsep//trim(EnvelopPrefix_post), & ! File containing all RCs corresponding to MCMC samples - teta uncertainty + remnant std uncertainty
                            SaveSpag=SaveSpag_post,&
                            SaveEnvelop=SaveEnvelop_post,&
                            propagationMatrix=PropMatrix_post,& ! let the user define what is propagated and what is not
                            HQFile=trim(workspace)//dirsep//trim(HQFile),& ! File containing H/Q obs and uncertainties/envelops
				            SaveHQ=SaveHQ,&
                            CookedMCMCFile=trim(workspace)//dirsep//trim(CookedFile),SaveCooked=SaveCooked,&
				            err=err,mess=mess)! error handling
    if(err>0) then; call BaRatin_ConsoleMessage(-6,trim(mess));endif
    call BaRatin_ConsoleMessage(7, '')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(DoH2Q) then
    call BaRatin_ConsoleMessage(8,'')
    call BaRatinConfig_Read_H2QPropagation(file=trim(workspace)//dirsep//trim(Config_H2Qpropagation),&
                                   HHfile=HHfile,&
                                   SaveSpag=SaveSpag_prop,SpagPrefix=SpagPrefix_prop,&
                                   SaveEnvelop=SaveEnvelop_prop,EnvelopPrefix=EnvelopPrefix_prop,&
                                   PropMatrix=PropMatrix_prop,&
                                   err=err,mess=mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
    call BaRatin_Propagation(MCMCFile=trim(workspace)//dirsep//trim(MCMCfile),&
                        Nburn=Nburn,Nread=Nread,Nslim=Nslim,&
                        Hfile=trim(HHFile),&
                        RC=RC,&
                        SpagFile=trim(workspace)//dirsep//trim(SpagPrefix_prop), &
                        EnvFile=trim(workspace)//dirsep//trim(EnvelopPrefix_prop), &
                        SaveSpag=SaveSpag_prop,&
                        SaveEnvelop=SaveEnvelop_prop,&
                        propagationMatrix=PropMatrix_prop,& ! let the user define what is propagated and what is not
				        err=err,mess=mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
    call BaRatin_ConsoleMessage(9, '')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Good bye!
call BaRatin_ConsoleMessage(999, '')
!read(*,*)

end program main

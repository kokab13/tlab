#include "dns_const.h"
#include "dns_error.h"

subroutine SCAL_READ_LOCAL(inifile)

    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS, only: efile, lfile, wfile, MAX_NSP
    use TLAB_VARS, only: inb_scal
    use TLAB_VARS, only: sbg
    use TLAB_PROCS
    use SCAL_LOCAL

    implicit none

    character*(*) inifile

! -------------------------------------------------------------------
    real(cp) dummy(MAX_NSP)
    integer(ci) idummy
    character*512 sRes
    character*32 bakfile

! ###################################################################
    bakfile = TRIM(ADJUSTL(inifile))//'.bak'

    call TLAB_WRITE_ASCII(lfile, 'Reading local input data')

! ###################################################################
! Scalar fluctuation field
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[IniFields]')
    call TLAB_WRITE_ASCII(bakfile, '#Scalar=<option>')
    call TLAB_WRITE_ASCII(bakfile, '#ThickIniS=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#ProfileIniS=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#YCoorIniS=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#NormalizeS=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#Mixture=<string>')

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'Scalar', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; flag_s = 0
    else if (TRIM(ADJUSTL(sRes)) == 'layerbroadband') then; flag_s = 1
    else if (TRIM(ADJUSTL(sRes)) == 'layerdiscrete') then; flag_s = 2
!  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'both'           ) THEN; flag_s = 3
    else if (TRIM(ADJUSTL(sRes)) == 'planebroadband') then; flag_s = 4
    else if (TRIM(ADJUSTL(sRes)) == 'planediscrete') then; flag_s = 5
    else if (TRIM(ADJUSTL(sRes)) == 'deltabroadband') then; flag_s = 6
    else if (TRIM(ADJUSTL(sRes)) == 'deltadiscrete') then; flag_s = 7
    else if (TRIM(ADJUSTL(sRes)) == 'fluxbroadband') then; flag_s = 8
    else if (TRIM(ADJUSTL(sRes)) == 'fluxdiscrete') then; flag_s = 9; end if

    Sini(:) = sbg(:) ! default geometry and scaling of perturbation

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'ProfileIniS', 'GaussianSurface', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; Sini(:)%type = PROFILE_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; Sini(:)%type = PROFILE_GAUSSIAN
    else if (TRIM(ADJUSTL(sRes)) == 'gaussiansurface') then; Sini(:)%type = PROFILE_GAUSSIAN_SURFACE
    else if (TRIM(ADJUSTL(sRes)) == 'gaussiansymmetric') then; Sini(:)%type = PROFILE_GAUSSIAN_SYM
    else if (TRIM(ADJUSTL(sRes)) == 'gaussianantisymmetric') then; Sini(:)%type = PROFILE_GAUSSIAN_ANTISYM
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Wrong ProfileIni parameter.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniS', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &    ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then; Sini(:)%thick = sbg(:)%thick; 
    else
        dummy = 0.0_cp; idummy = MAX_NSP
        call LIST_REAL(sRes, idummy, dummy)
        Sini(:)%thick = dummy(:)
        if (idummy /= inb_scal) then         ! Consistency check
            if (idummy == 1) then
                Sini(2:)%thick = Sini(1)%thick
                call TLAB_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using ThickIniS(1) for all scalars.')
            else
                call TLAB_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. ThickIniS size does not match number of scalars.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if
    end if

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIniS', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &    ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then; Sini(:)%ymean = sbg(:)%ymean; 
    else
        dummy = 0.0_cp; idummy = MAX_NSP
        call LIST_REAL(sRes, idummy, dummy)
        Sini(:)%ymean = dummy(:)
        if (idummy /= inb_scal) then         ! Consistency check
            if (idummy == 1) then
                Sini(2:)%ymean = Sini(1)%ymean
                call TLAB_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using YCoorIniS(1) for all scalars.')
            else
                call TLAB_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. YCoorIniS size does not match number of scalars.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if
    end if

    ! do idummy = 1, inb_scal; Sini(idummy)%parameters(:) = 0.0_cp; end do

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'NormalizeS', '-1.0', sRes)
    norm_ini_s(:) = 0.0_cp; idummy = MAX_NSP
    call LIST_REAL(sRes, idummy, norm_ini_s)
    if (idummy /= inb_scal) then            ! Consistency check
        if (idummy == 1) then
            norm_ini_s(2:) = norm_ini_s(1)
            call TLAB_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using NormalizeS(1) for all scalars.')
        else
            call TLAB_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. NormalizeS size does not match number of scalars.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if
    end if

    call SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeR', '0.0', norm_ini_radiation) ! Radiation field

! Additional parameters
    call SCANINICHAR(bakfile, inifile, 'IniFields', 'Mixture', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; flag_mixture = 0
    else if (TRIM(ADJUSTL(sRes)) == 'equilibrium') then; flag_mixture = 1
    else if (TRIM(ADJUSTL(sRes)) == 'loadfields') then; flag_mixture = 2
    end if

! ###################################################################
! Discrete Forcing
! ###################################################################
    call DISCRETE_READBLOCK(bakfile, inifile, 'Discrete', fp)
    ! Modulation type in fp%type

!   specific for this tool
    call SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))
    call SCANINIREAL(bakfile, inifile, 'Discrete', 'ThickStep', '-1.0', fp%parameters(2))

    return
end subroutine SCAL_READ_LOCAL

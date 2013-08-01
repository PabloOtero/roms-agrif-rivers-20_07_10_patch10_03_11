!
! $Id: get_river.F,v 1.0 2010/05/04 13:56:02 potero $
!
#include "cppdefs.h"

#if defined PSOURCE
                              ! Read river data
      subroutine get_river   ! at the appropriate
                              ! time from forcing netCDF file. 

      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "forces.h"
# include "netcdf.inc"
# include "sources.h"

#if defined ANA_PSOURCE
      real cff, buff(Msrc)
      integer is, i, j, ierr, lstr, lvar,
     &    lenstr, nf_fread, advance_cycle
      integer startriv(2), countriv(2)
!
! Initialization: Inquire about the contents of forcing netCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      if (itriv.eq.0 .or. iic.eq.0) then
        lstr=lenstr(riverfile)
c**        call opencdf (riverfile,N)
        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing netCDF file for reading. Find and
! save IDs for relevant variables, determine whether river data
! components exist as matrix (varndims=2), vectors (1) or scalars (0).
!
        if (ncidriv.eq.-1) then
          ierr=nf_open (riverfile(1:lstr), nf_nowrite, ncidriv)
          if (ierr .ne. nf_noerr) goto 4                 !--> ERROR 
        endif

        ierr=nf_inq_varid (ncidriv, 'river_time', riv_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'river_time', riverfile(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxRUNOFF))
        ierr=nf_inq_varid (ncidriv,vname(1,indxRUNOFF)(1:lvar),runoffid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidriv, runoffid, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lsrivgrd=1
            else
              lsrivgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxRUNOFF)(1:lvar),riverfile(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxTRIVER))
        ierr=nf_inq_varid (ncidriv,vname(1,indxTRIVER)(1:lvar),triverid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidriv, triverid, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lsrivgrd=1
            else
              lsrivgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxTRIVER)(1:lvar),riverfile(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxSRIVER))
        ierr=nf_inq_varid (ncidriv,vname(1,indxSRIVER)(1:lvar),sriverid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidriv, sriverid, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lsrivgrd=1
            else
              lsrivgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxSRIVER)(1:lvar),riverfile(1:lstr)
          goto 99                                         !--> ERROR
        endif

! Determine whether there is cycling to reuse the input data and
! what is cycling period "riv_cycle", find initial cycling index
! "riv_ncycle" and record index "riv_rec". Set initial value for
! time index "itriv" and both time record bounds to large negative
! artificial values, so that it will trigger the logic in reading
! part below.
!
        call set_cycle (ncidriv,   riv_tid,    ntriv,
     &                  riv_cycle, riv_ncycle, riv_rec)

        if (may_day_flag.ne.0) return    !-->  EXIT
        itriv=2
        riv_time(1)=-1.E+20
        riv_time(2)=-1.E+20
      endif                   !<-- itsms.eq.0 .or. iic.eq.0

!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cycling indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!

  1   i=3-itriv
       cff=time+0.5*dt  

        if (riv_time(i).le.cff .and. cff.lt.riv_time(itriv)) return

        ierr=advance_cycle (riv_cycle, ntriv, riv_ncycle, riv_rec)
        if (ierr.ne.0) goto 7                           !--> ERROR

        ierr=nf_get_var1_FTYPE(ncidriv, riv_tid, riv_rec, cff)
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'riv_time', riv_rec 
          goto 99                                         !--> ERROR
        endif

        riv_time(i)=cff*day2sec+riv_cycle*riv_ncycle
        if (riv_time(itriv).eq.-1.E+20) riv_time(itriv)=riv_time(i)
        if (lsrivgrd.eq.1) then
          ierr=nf_inq_dimid (ncidriv, 'river_id', j)
          ierr=nf_inq_dimlen(ncidriv,j,j)    
          startriv(1)=1
	  startriv(2)=riv_rec
	  countriv(1)=j
          countriv(2)=1

          ierr=nf_get_vara_double(ncidriv,runoffid,startriv,countriv,buff)
          runofftrp(1:Msrc,i)=buff
          ierr=nf_get_vara_double(ncidriv,triverid,startriv,countriv,buff)
          trivertrp(1:Msrc,i)=buff
          ierr=nf_get_vara_double(ncidriv,sriverid,startriv,countriv,buff)
          !srivertrp(1:Msrc,i)=buff
          srivertrp(1:Msrc,i)=5.
        else
          ierr=nf_get_var1_FTYPE(ncidriv,runoffid,riv_rec,runofftrp(1:Msrc,i))
          ierr=nf_get_var1_FTYPE(ncidriv,triverid,riv_rec,trivertrp(1:Msrc,i))
          ierr=nf_get_var1_FTYPE(ncidriv,sriverid,riv_rec,srivertrp(1:Msrc,i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'runofftr', riv_rec 
          goto 99                                         !--> ERROR
        endif

       itriv=i

       write(stdout,'(6x,A,1x,A,1x,g12.5,1x,I4)') 'GET_RIVER --',
     &              'Read river discharge data for time =', cff
       PRINT*
       PRINT*, '--DISCHARGE (m3/s)----Temperature(ºC)----Salinity--'
       do is=1,Nsrc
       write(stdout,'(1x,A,4x,f6.2,4x,f6.2,4x,f6.2,4x,I4)') nameriver(is), 
     &         runofftrp(Ident(is),i), trivertrp(Ident(is),i),
     &         srivertrp(Ident(is),i)
       enddo
       PRINT*

#ifdef MPI
     &                                                      , mynode
#endif
      if (ntriv.gt.1) goto 1
      if (ntriv.eq.1) return
!
! Sort out error messages: The following portion of the code is 
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_RIVER - unable to find river variable: ',A,
     &                          /,14x,'in forcing netCDF file: ',A)
      goto 99
  4   write(stdout,5) riverfile(1:lstr)
  5   format(/,' GET_RIVER - unable to open river netCDF file:',
     &                                                        1x,A)
  6   format(/,' GET_RIVER - error while reading variable: ',A,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) riv_rec, ntriv, riverfile(1:lstr), tdays,
     &                                riv_time(itriv)*sec2day
  8   format(/,' GET_RIVER - ERROR: requested time record ',I4,
     &       1x,'exceeds the last available', /,14x,'record ',I4,
     &       1x,'in river netCDF file: ',a, /,14x,'TDAYS = ',
     &       g12.4,2x,'last available RIVER_TIME = ',g12.4)
 99   may_day_flag=2
      return 
#endif /* !ANA_PSOURCE*/
      end
#endif /* !PSOURCE*/

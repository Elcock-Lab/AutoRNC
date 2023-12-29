
C this code reads a NASCENT_CTERMINUS.pdb and a RIBOSOME.pdb
C and adds new blocking HETATM entries in a separate file that
C can later be catted to the end of the RIBOSOME.pdb

C blocking atoms are placed around the final residue of the
C NASCENT_CTERMINUS.pdb with all points that are closer to the N atom
C than the C atom removed - this should allow the nascent chain to build
C forwards...

      program add_some_blocking_atoms

      integer iatom,ile,ilr,ires,ilast,ihave
      character linia*76,head*30,tail*22,subunit*1,junk*80
      character string*4,atom_n*5,res_n*4,strg*4
      character*80 filein_CTER
      character*80 filein_RIBO
      character*80 fileout
      character*66 char66
      character*4  char4

      real x_CTER(1000000)
      real y_CTER(1000000)
      real z_CTER(1000000)
      real x_RIBO(1000000)
      real y_RIBO(1000000)
      real z_RIBO(1000000)

      real*4 x,y,z

      call getarg(1,junk)
      read(junk,*)filein_CTER
      call getarg(2,junk)
      read(junk,*)filein_RIBO
      call getarg(3,junk)
      read(junk,*)fileout
      call getarg(4,junk)
      read(junk,*)radius ! radius around C atom
      call getarg(5,junk)
      read(junk,*)clash ! omit points closer to RIBO than this
      call getarg(6,junk)
      read(junk,*)iblock! which res of Nascent chain to do blocking at

      radius2=radius*radius
      clash2=clash*clash
 
      natm_CTER=0
      nblockN=0
      nblockC=0

      x_N=-999.9
      y_N=-999.9
      z_N=-999.9
      x_C=-999.9
      y_C=-999.9
      z_C=-999.9
      open(11,file=filein_CTER,status='old')

10    read(11,'(a66)',end=20)char66
      if(char66(1:4).eq.'ATOM'.or.char66(1:4).eq.'HETA') then
        natm_CTER=natm_CTER+1
        read(char66,'(12x,a4,14x,3f8.3)')char4,xt,yt,zt
        x_CTER(natm_CTER)=xt 
        y_CTER(natm_CTER)=yt 
        z_CTER(natm_CTER)=zt 

C we check if this is a N or C - every time we find one we increment the
C counters and we keep them if we're still looking for the iblock'th
C residue - note we use ".le." here instead of ".eq." to handle cases
C where a user doesn't know what residue to use and so just puts in a
C big number - the code as written will auto-keep the last one

        if(char4.eq.' N  ') then
          nblockN=nblockN+1
          if(nblockN.le.iblock) then
            x_N=xt
            y_N=yt
            z_N=zt
          endif
        elseif(char4.eq.' C  ') then
          nblockC=nblockC+1
          if(nblockC.le.iblock) then
            x_C=xt
            y_C=yt
            z_C=zt
          endif
        endif

      endif
      goto 10
20    close(11)

C as a sanity check, write out the coords of the N and C atoms

      write(*,*)'x_N etc ',x_N,y_N,z_N
      write(*,*)'x_C etc ',x_C,y_C,z_C

C now read the RIBOSOME.pdb

      natm_RIBO=0
      open(11,file=filein_RIBO,status='old')

30    read(11,'(a66)',end=40)char66
      if(char66(1:4).eq.'ATOM'.or.char66(1:4).eq.'HETA') then
        natm_RIBO=natm_RIBO+1
        read(char66,'(12x,a4,14x,3f8.3)')char4,xt,yt,zt
        x_RIBO(natm_RIBO)=xt
        y_RIBO(natm_RIBO)=yt
        z_RIBO(natm_RIBO)=zt
      endif
      goto 30
40    close(11)

      write(*,*)'natm_CTER= ',natm_CTER
      write(*,*)'natm_RIBO= ',natm_RIBO

C now get ready to write out 

      open(11,file=fileout,status='unknown')

C loop in x,y,z - assume 1A spacing

      nout=0

      ilim=int(radius)+1

      do i=-ilim,+ilim
        do j=-ilim,+ilim
          do k=-ilim,+ilim

C skip this point if radius exceeds input value

            rad_tmp2=real(i)**2+real(j)**2+real(k)**2
            if(rad_tmp2.gt.radius2) cycle

            xp=x_C+real(i)
            yp=y_C+real(j)
            zp=z_C+real(k)

C skip this point if closer to CTER N than CTER C

            dist2N=(x_N-xp)**2+
     &             (y_N-yp)**2+
     &             (z_N-zp)**2
            dist2C=(x_C-xp)**2+
     &             (y_C-yp)**2+
     &             (z_C-zp)**2
            if(dist2N.le.dist2C) cycle

C assume for now that we will keep this grid point

            ikeep=1

C skip this point if it clashes with CTER

            do n=1,natm_CTER
              dist2=((x_CTER(n)-xp)**2+
     &               (y_CTER(n)-yp)**2+
     &               (z_CTER(n)-zp)**2)
              if(dist2.le.clash2) then
                ikeep=0
                exit
              endif
            enddo

            if(ikeep.eq.0) cycle

C skip this point if it clashes with RIBO

            do n=1,natm_RIBO
              dist2=((x_RIBO(n)-xp)**2+
     &               (y_RIBO(n)-yp)**2+
     &               (z_RIBO(n)-zp)**2)
              if(dist2.le.clash2) then
                ikeep=0
                exit
              endif
            enddo

            if(ikeep.eq.0) cycle

C otherwise, write it out...

            nout=nout+1
            nout2=mod(nout,10000)
            write(11,50)nout,nout2,xp,yp,zp
50    format('HETATM',i5,'  N   XXX  ',i4,4x,3f8.3,'  1.00 -1.00')

          enddo
        enddo
      enddo

      close(11)

      write(*,*)'wrote #blocking atoms= ',nout

      stop
      end


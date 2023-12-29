
C this code is for use with AutoRNC - it'll read in a protein pdb and
C depending on what setting for sidechain_mode - if 1 it'll keep only 
C the CB of the sidechain - if 2 it'll keep no part of the protein 
C sidechain

      program transform

      integer iatom,ile,ilr,ires,ilast,ihave
      character linia*76,head*30,tail*22,subunit*1,junk*80
      character string*4,atom_n*5,res_n*4,strg*4
      character*80 filein_CTER
      character*80 filein
      character*80 fileout
      character*66 char66
      character*4  char4

      real x_CTER(1000000)
      real y_CTER(1000000)
      real z_CTER(1000000)
      real x_RIBO(1000000)
      real y_RIBO(1000000)
      real z_RIBO(1000000)

      integer sidechain_mode

      real*4 x,y,z

      call getarg(1,junk)
      read(junk,*)filein
      call getarg(2,junk)
      read(junk,*)fileout
      call getarg(3,junk)
      read(junk,*)sidechain_mode

      natm_RIBO=0
      open(11,file=filein,status='unknown')
      open(12,file=fileout,status='unknown')

30    read(11,'(a66)',end=40)char66
      if(char66(1:4).eq.'ATOM'.or.char66(1:4).eq.'HETA') then
        if(char66(18:20).eq.'ALA'.or.
     &     char66(18:20).eq.'ARG'.or.
     &     char66(18:20).eq.'ASN'.or.
     &     char66(18:20).eq.'ASP'.or.
     &     char66(18:20).eq.'CYS'.or.
     &     char66(18:20).eq.'GLN'.or.
     &     char66(18:20).eq.'GLU'.or.
     &     char66(18:20).eq.'GLY'.or.
     &     char66(18:20).eq.'HIS'.or.
     &     char66(18:20).eq.'ILE'.or.
     &     char66(18:20).eq.'LEU'.or.
     &     char66(18:20).eq.'LYS'.or.
     &     char66(18:20).eq.'MET'.or.
     &     char66(18:20).eq.'MSE'.or.
     &     char66(18:20).eq.'PHE'.or.
     &     char66(18:20).eq.'PRO'.or.
     &     char66(18:20).eq.'SER'.or.
     &     char66(18:20).eq.'THR'.or.
     &     char66(18:20).eq.'TRP'.or.
     &     char66(18:20).eq.'TYR'.or.
     &     char66(18:20).eq.'VAL') then
          if(sidechain_mode.eq.1) then
            if(char66(13:16).ne.' N  '.and.
     &         char66(13:16).ne.' CA '.and.
     &         char66(13:16).ne.' CB '.and.
     &         char66(13:16).ne.' C  '.and.
     &         char66(13:16).ne.' O  ') then
            else
              write(12,'(a66)')char66
            endif
          elseif(sidechain_mode.eq.2) then
            if(char66(13:16).ne.' N  '.and.
     &         char66(13:16).ne.' CA '.and.
     &         char66(13:16).ne.' C  '.and.
     &         char66(13:16).ne.' O  ') then
            else
              write(12,'(a66)')char66
            endif
          endif
        else
          write(12,'(a66)')char66
        endif
      else
        write(12,'(a66)')char66
      endif
      goto 30
40    close(11)
      close(12)

      stop
      end


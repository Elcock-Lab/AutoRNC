
C this code is for use with AutoRNC - it'll read in a pdb file - e.g.
C one from AlphaFold2 and remove any hydrogens, as defined by entries
C that have a "H" in column 78...

      program transform

      integer iatom,ile,ilr,ires,ilast,ihave
      character linia*76,head*30,tail*22,subunit*1,junk*80
      character string*4,atom_n*5,res_n*4,strg*4
      character*80 filein_CTER
      character*80 filein
      character*80 fileout
      character*80 char80
      character*4  char4

      real x_CTER(1000000)
      real y_CTER(1000000)
      real z_CTER(1000000)
      real x_RIBO(1000000)
      real y_RIBO(1000000)
      real z_RIBO(1000000)

      real*4 x,y,z

      call getarg(1,junk)
      read(junk,*)filein
      call getarg(2,junk)
      read(junk,*)fileout

      natm_RIBO=0
      open(11,file=filein,status='unknown')
      open(12,file=fileout,status='unknown')

30    read(11,'(a80)',end=40)char80
      if(char80(78:78).eq.'H') goto 30
      write(12,'(a80)')char80
      goto 30
40    close(11)
      close(12)

      stop
      end


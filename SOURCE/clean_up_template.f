
C this is simple code that overwrites the occupancy and beta factors for
C ATOM entries in a .pdb file - this is a preprocessing step in the code

      program clean_template
 
      implicit double precision(a-h,o-z)
      character*80 fname1
      character*80 fname2
      character*80 junk
      character*19 uplbl
      character*10 nxtlbl
      character*80 string
      character*3  atmname
      call getarg(1,junk)
      read(junk,*)fname1
      call getarg(2,junk)
      read(junk,*)fname2

      open(unit=10,file=fname1,status="unknown")
      open(unit=11,file=fname2,status="unknown")

16    read(10,20,end=30)string
      if(string(1:4).ne.'ATOM') then
        write(11,20)string
      else
        string(55:66)='  1.00 -1.00'
        write(11,'(a66)')string(1:66)
      endif

      goto 16

20    format(a80)
22    format(a22,i4,a40)
30    close(10)
      close(11)
      stop
      end

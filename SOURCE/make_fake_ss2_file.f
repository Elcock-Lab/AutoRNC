
C this code reads a list of secondary structure elements and a fasta and
C writes out a fake .ss2 file that is used by ahemodel_2023_v1.2.exe
C later to force the construction of helices and strands at positions
C specified by the user...

      program make_fake_ss2
 
      implicit real(a-h,o-z)
      parameter   (pi=3.141592654)
      character*1 c          
      character*3 a        
      character*3 r        
      integer     n,num_box
      real, allocatable :: xs(:)
      real, allocatable :: ys(:)
      real, allocatable :: zs(:)

      type gridtype
        integer                     :: num
        real, dimension(:), pointer :: x
        real, dimension(:), pointer :: y
        real, dimension(:), pointer :: z
        integer, dimension(:), pointer :: i ! atom number
      end type gridtype
      type (gridtype), allocatable  :: grid (:,:,:)    
      type (gridtype)  :: temp_arr                    

      character*80 junk        

      character*1  char1
      character*60 char60
      character*80 char80
      character*80 list_of_ss_elems
      character*80 fasta_file
      character*80 ss2_file
      character*1 seq(100000)
      integer ibeg(100000)
      integer iend(100000)
      integer ityp(100000)

      call getarg(1,junk)
      read(junk,*)list_of_ss_elems   ! input file of SS elems
      call getarg(2,junk)
      read(junk,*)fasta_file         ! input fasta file
      call getarg(3,junk)
      read(junk,*)nres               ! input #res
      call getarg(4,junk)
      read(junk,*)ss2_file           ! output ss2 file                    

      do i=1,nres
        ityp(i)=0 ! assume all res coil to start
      enddo

C read list of secondary structure elements

      open(unit=11,file=list_of_ss_elems,status='unknown')
      read(11,*)
10    read(11,'(a80)',end=20)char80
      read(char80,*)l,char1,i,j
      write(*,*)'l etc ',l,i,j,char1
      if(char1.eq.'H') then
        ityp(i:j)=1 
      elseif(char1.eq.'E') then
        ityp(i:j)=2
      endif
      goto 10
20    close(11) 

C read fasta file - assume blocks of 60 characters...

      mres=0
      open(unit=11,file=fasta_file,status='unknown')
      read(11,*)
30    read(11,'(a60)',end=40)char60
      do i=1,60
        if(char60(i:i).eq.' ') goto 40
        mres=mres+1
        seq(mres)=char60(i:i)
      enddo
      goto 30
40    close(11)
 
      open(unit=11,file=ss2_file,status='unknown')
      write(11,'("# PSIPRED VFORMAT (PSIPRED V4.0)")')
      write(11,'(1x)')
      do i=1,nres
        f1=0.0
        f2=0.0
        f3=0.0
        if(ityp(i).eq.0) then
          char1='C'
          f1=1.0
        elseif(ityp(i).eq.1) then
          char1='H'
          f2=1.0
        elseif(ityp(i).eq.2) then
          char1='E'
          f3=1.0
        endif
        write(11,50)i,seq(i),char1,f1,f2,f3
50      format(i4,1x,a,1x,a,1x,3f7.3)
      enddo
      close(11)
 
      stop
      end

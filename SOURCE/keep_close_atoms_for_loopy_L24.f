
C this code reads in a L24 pdb and a complete ribosome pdb and writes
C out only those ribosome atoms that are close to the L24 in a new pdb
C file - note that it gives them fake GLY names to make sure that loopy
C behaves okay later...

      program do_some_stuff_with_l24
 
      implicit real (a-h,o-z)
      common/aherandom/idumahe
      integer nmol1(100000)
      integer nmol2(100000)
      integer iatm1(100000)
      integer iatm2(100000)
      real x1(100000)
      real y1(100000)
      real z1(100000)
      real x2(300000)
      real y2(300000)
      real z2(300000)
      character*3 jnk,knk,nam        
      character*1 lnk,lnk_lst
      character*50 junk
      character*50 file1,file1pseudo
      character*50 file2,file2pseudo
      character*50 filefinal
      character*80 string   
      character*80 strong(300000)
      character*1  chainID
      character*4  char4
      character*6  char6
      character*1  char1
      parameter                                                     
     &         (pi=3.141592654,                                        
     &       twopi=6.283185307,                                       
     &         eps=0.001,                                               
     &        eps2=0.999999,                                           
     &        eps3=0.0001,                                             
     &      dtorad=0.0174532925)   
      call getarg(1,junk)
      read(junk,*) file1    ! L24 pdb file
      call getarg(2,junk)
      read(junk,*) file2    ! ribosome pdb file
      call getarg(3,junk)
      read(junk,*) filefinal! output pdb file
      call getarg(4,junk)
      read(junk,*) cut_dist ! distance cutoff
      call getarg(5,junk)
      read(junk,*) char1    ! chainID for fake gly
  

      cut_dist2=cut_dist**2
      natm1=0
      natm2=0
      open(unit=10,file=file1,status='unknown')
      open(unit=11,file=file2,status='unknown')
      open(unit=12,file=filefinal,status='unknown')
10    read(10,12,end=15)string
12    format(a80)
      if(string(1:4).ne.'ATOM') goto 10
      natm1=natm1+1
      read(string,14)x,y,z
      read(string,'(21x,a1)')chainID
14    format(30x,3f8.3)
      x1(natm1)=x
      y1(natm1)=y
      z1(natm1)=z
      goto 10
15    close(10)
20    read(11,12,end=25)string
      if(string(1:4).ne.'ATOM'.and.string(1:4).ne.'HETA') goto 20

C skip those atoms that are part of L24 already

      if(string(22:22).eq.chainID) goto 20
      natm2=natm2+1
      read(string,14)x,y,z
      x2(natm2)=x
      y2(natm2)=y
      z2(natm2)=z
      strong(natm2)=string
      goto 20
25    close(11)

      write(*,*)'natm1 = ',natm1
      write(*,*)'natm2 = ',natm2

      ngly=1
      ningly=0
      natm=0
      do j=1,natm2
        do i=1,natm1
          dist2=(x1(i)-x2(j))**2+
     &          (y1(i)-y2(j))**2+
     &          (z1(i)-z2(j))**2
          if(dist2.le.cut_dist2) then
            ningly=ningly+1
            if(ningly.gt.4) then
              ngly=ngly+1
              ningly=1
            endif
            write(char4,'(i4)')ngly
            strong(j)(23:26)=char4
            if(ningly.eq.1) then
              strong(j)(13:16)=' N  '
            elseif(ningly.eq.2) then
              strong(j)(13:16)=' CA '
            elseif(ningly.eq.3) then
              strong(j)(13:16)=' C  '
            elseif(ningly.eq.4) then
              strong(j)(13:16)=' O  '
            endif
            strong(j)(22:22)=char1
            strong(j)(1:6)='ATOM  '
            strong(j)(18:20)='GLY'
            natm=natm+1
            write(char6,'(i6)')natm
            strong(j)(6:11)=char6
            write(*,*)'wrote atom # ',j,i,sqrt(dist2)
            write(12,'(a80)')strong(j)(1:80)
            jlast=j
            exit
          endif
        enddo
      enddo

C need to remember to write out a complete GLY residue otherwise loopy
C throws a fit if I remember correctly...

      ningly_final=ningly
      do ningly=ningly_final+1,4
        if(ningly.eq.1) then
          strong(jlast)(13:16)=' N  '
        elseif(ningly.eq.2) then
          strong(jlast)(13:16)=' CA '
        elseif(ningly.eq.3) then
          strong(jlast)(13:16)=' C  '
        elseif(ningly.eq.4) then
          strong(jlast)(13:16)=' O  '
        endif
        strong(jlast)(22:22)=char1
        strong(jlast)(1:6)='ATOM  '
        strong(jlast)(18:20)='GLY'
        natm=natm+1
        write(char6,'(i6)')natm
        strong(jlast)(6:11)=char6
        write(12,'(a80)')strong(jlast)(1:80)
      enddo

      stop
      end

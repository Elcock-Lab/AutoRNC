
C this is a preprocessing code that reads in a full-length protein .pdb
C file and keeps only those residues specified by the user

      program cull_protein
 
      implicit double precision(a-h,o-z)
      character*80 fname1
      character*80 fname2
      character*80 fname3
      character*50 junk
      character*19 uplbl
      character*10 nxtlbl
      character*80 char80
      character*10000 char10000
      character*3  char3
      character*5  char5
      character*94 char94
      integer      iwant(10000)
      character*28 my_fmt

      call getarg(1,junk)
      read(junk,*)fname1 ! input pdb file
      call getarg(2,junk)
      read(junk,*)fname3 ! output pdb file
      call getarg(3,junk)
      read(junk,*)ibeg   ! first residue to keep
      call getarg(4,junk)
      read(junk,*)iend   ! last residue to keep
      call getarg(5,junk)
      read(junk,*)nfinal ! final residue of construct

      open(unit=11,file=fname1,status="unknown")

C first figure out how many residues there are in the protein

      nres=0
10    read(11,'(a80)',end=20)char80
      if(char80(1:4).ne.'ATOM') goto 10
      read(char80,'(22x,i4)')ires
      nres=max(nres,ires)
      goto 10
20    rewind(11)

      write(*,*)"# of total residues = ",nres

C now set iwant=1 for those residues from ibeg to iend

      iwant=0
      iwant(ibeg:iend)=1
      nres_actual=nfinal

C now go and extract these residues from the input pdb and put them in
C the output pdb - note that we write out all non-ATOM entries
C automatically...

      nseqres=0
      idoneseqres=0
      open(unit=12,file=fname3,status="unknown")
30    read(11,'(a80)',end=40)char80
      if(char80(1:4).ne.'ATOM') then
        if(char80(1:6).eq.'SEQRES') then
          nseqres=nseqres+1
          if(nseqres*13.le.nres_actual) then
            write(12,'(a80)')char80
          elseif(idoneseqres.eq.0) then
            ntodo=nres_actual-(nseqres-1)*13
            nbeg=18+4*ntodo+1
            nend=70
            do n=nbeg,nend
              char80(n:n)=' '
            enddo
            write(12,'(a80)')char80
            idoneseqres=1
          endif
          goto 30
        else     
          write(12,'(a80)')char80
          goto 30
        endif
      endif

C if we get here then it's an atom entry...

      read(char80,'(22x,i4)')ires
      if(iwant(ires).eq.0) goto 30
      write(12,'(a80)')char80
      goto 30
40    rewind(11)

      stop
      end

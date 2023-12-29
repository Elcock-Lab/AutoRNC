
C this code performs a preprocessing step for AutoRNC - it generates the
C TEMPLATE .pdb file and a corresponding alignment file that are read by
C ahemodel when it begins the chain construction procedure...

      program prepare_template_for_ahemodel
 
      integer iatom,ile,ilr,ires,ilast,ihave,num_trans
      integer num_mons_orig
      character*12 my_fmt
      character*1000 dipeptide_path
      character*1000 dipepname
      dimension x1(100),x2(100),x3(100),x4(100)
      dimension y1(100),y2(100),y3(100),y4(100)
      dimension z1(100),z2(100),z3(100),z4(100)
      character linia*76,head*30,tail*22,subunit*1,junk*22
      character string*4,atom_n*5,res_n*4,strg*4
      character*2000 fasta
      character*52 chain   
      character*1  chainID
c     character*60 seq_chunk
      character*10000 seq_chunk
      character*10000 seq_of_fasta
      character*10000 seq2
      character*10000 seq3
      character*70 outstring
      character*3  rnam
      character*3  rnam_final(1000000)
      character*4  char4
      character*80 crap,fasta_in,alignment_out,template_out,template_in
      character*66 string66

      real   x,y,z
      call getarg(1,crap)
      read(crap,*)fasta_in 
      call getarg(2,crap)
      read(crap,*)template_in ! NASCENT
      call getarg(3,crap)
      read(crap,*)alignment_out
      call getarg(4,crap)
      read(crap,*)template_out
      call getarg(5,crap)
      read(crap,*)nres_nascent ! NASCENT
      call getarg(6,crap)
      read(crap,*)nstart       ! 1st res to build in N-ter direction

C NASCENT code: set iter=2 so we build from the C-terminus forwards

      iter=2

      call init_random_seed()

C here is how we *should* be able to define the path...
C it doesn't work so we hard code it below...
 
c     dipeptide_path=adjustl(dipeptide_path)
c     klm=len(trim(dipeptide_path))
c     write(*,*)'klm = ',klm
c     write(*,*)'dipeptide_path = ',dipeptide_path(1:20)

      dipeptide_path='/home/aelcock/2017_ECOLI_CELL/PISCES'
      klm=36

      do i=1,10000
        seq_of_fasta(i:i)=' '
        seq2(i:i)=' '
        seq3(i:i)=' '
      enddo

199   format(a80)
c
      open(2,file=fasta_in,status='old')
      open(3,file=template_out,status='unknown')
      open(4,file=alignment_out,status='unknown')

      if(iter.eq.1) then
        write(3,'("REMARK_AHE0: FAKE UNFOLDED N-ter TEMPLATE")')
      else
        write(3,'("REMARK_AHE0: FAKE UNFOLDED C-ter TEMPLATE")')
      endif

C read the fasta files for the input sequences

      nres=0
      read(2,*) ! skip top line
26    read(2,27,end=28)seq_chunk
27    format(a10000)

      do i=1,len(trim(seq_chunk)) ! line don't need to have 60 chars
        if(seq_chunk(i:i).ne.' ') then
          nres=nres+1
          seq_of_fasta(nres:nres)=seq_chunk(i:i)
        else
          goto 28
        endif
      enddo

      goto 26
28    close(2)

      chainID='A'

      write(*,*)'nres = ',nres
      do i=1,70
        outstring(i:i)=' '
      enddo
      outstring(1:19)='SEQRES   1 X  123  '
      outstring(12:12)=chainID
ccc   write(outstring(13:17),'(i5)')nres
      write(outstring(13:17),'(i5)')1
      nline=0

      if(iter.eq.1) then
        nbeg=1
        nend=1
      elseif(iter.eq.2) then
        nbeg=nres_nascent ! NASCENT
        nbeg=1
        nend=nres_nascent 
      endif

ccc   do mp=1,nres,13
      do mp=nbeg,nend,13
        nline=nline+1
        write(outstring(7:10),'(i4)')nline
        ibeg=mp
ccc     iend=min(nres,mp+12)
        iend=min(nend,mp+12)
        j=0
        write(*,*)'ibeg,iend ',ibeg,iend
        do i=ibeg,iend
          j=j+1
          if(seq_of_fasta(i:i).eq.'A') rnam='ALA'
          if(seq_of_fasta(i:i).eq.'R') rnam='ARG'
          if(seq_of_fasta(i:i).eq.'N') rnam='ASN'
          if(seq_of_fasta(i:i).eq.'D') rnam='ASP'
          if(seq_of_fasta(i:i).eq.'C') rnam='CYS'
          if(seq_of_fasta(i:i).eq.'Q') rnam='GLN'
          if(seq_of_fasta(i:i).eq.'E') rnam='GLU'
          if(seq_of_fasta(i:i).eq.'G') rnam='GLY'
          if(seq_of_fasta(i:i).eq.'H') rnam='HIS'
          if(seq_of_fasta(i:i).eq.'I') rnam='ILE'
          if(seq_of_fasta(i:i).eq.'L') rnam='LEU'
          if(seq_of_fasta(i:i).eq.'K') rnam='LYS'
          if(seq_of_fasta(i:i).eq.'M') rnam='MET'
          if(seq_of_fasta(i:i).eq.'F') rnam='PHE'
          if(seq_of_fasta(i:i).eq.'P') rnam='PRO'
          if(seq_of_fasta(i:i).eq.'S') rnam='SER'
          if(seq_of_fasta(i:i).eq.'T') rnam='THR'
          if(seq_of_fasta(i:i).eq.'W') rnam='TRP'
          if(seq_of_fasta(i:i).eq.'Y') rnam='TYR'
          if(seq_of_fasta(i:i).eq.'V') rnam='VAL'
          i1=20+(j-1)*4
          i2=22+(j-1)*4
          outstring(i1:i2)=rnam
          rnam_final(i)=rnam
cc        write(*,*)'i1,i2 ',i1,i2
        enddo

C remember to blank out the rest of the last line

        do j=i2+1,70
          outstring(j:j)=' '
        enddo

        write(3,'(a70)')outstring(1:70)
      enddo

C at this point we need to get coordinates for the C-ter most residue
C we just read these from template_in but we change the resname to that
C of the last residue desired from the mRNA (i.e. nres_nascent)
C for the "folded" code we also change the chainID to "A" and the res# to
C nres_nascent

      open(unit=11,file=template_in,status='unknown')

C we read through once just to find the highest res#
      ires_max=-999
4991  read(11,'(a66)',end=4992)string66
      if(string66(1:4).eq.'ATOM') then
        read(string66(23:26),'(i4)')ires
        ires_max=max(ires,ires_max)
      endif
      goto 4991
4992  rewind(11)

C and figure out the offset required to make each res have final#

      ires_mod=nres_nascent-ires_max

      nmods=0
      natms=0
5991  read(11,'(a66)',end=5992)string66
      if(string66(1:4).eq.'ATOM') then
        natms=natms+1
        string66(1:7)='ATOM   '
        write(char4,'(i4)')natms
        string66(8:11)=char4
        string66(22:22)='A'
        read(string66(23:26),'(i4)')ires
        write(char4,'(i4)')ires+ires_mod
        string66(23:26)=char4
        string66(18:20)=rnam_final(ires+ires_mod)

C skip writing this atom out if it's part of a residue that we want to
C rebuild - this helps us in cases where the construct's stall signal is
C shorter than the stall signal used in the ribosome crystal structure

        if(ires+ires_mod.le.nstart) goto 5991

      endif
      write(3,'(a66)')string66
      goto 5991
5992  close(11)
      close(3)

C now write out a completely fake alignment file for use in AHEMODEL

      write(4,'("########################################")')
      write(4,'("# Program: needle")')
      write(4,'("# Rundate: Wed  4 Sep 2019 13:54:09")')
      write(4,'("# Commandline: needle")')
      write(4,'("#    [-asequence] P0A6P9.fasta")')
      write(4,'("#    [-bsequence] TEMPLATE_ahemodel.001.fasta")')
      write(4,'("#    -gapopen 10.0")')
      write(4,'("#    -gapextend 0.5")')
      write(4,'("#    [-outfile] P0A6P9_TEMPLATE_ahemodel.001.out")')
      write(4,'("# Align_format: srspair")')
      write(4,'("# Report_file: P0A6P9_TEMPLATE_ahemodel.001.out")')
      write(4,'("########################################")')
      write(4,'("")')
      write(4,'("#=======================================")')
      write(4,'("#")')
      write(4,'("# Aligned_sequences: 2")')
      write(4,'("# 1: ENO_ECOLI")')
      write(4,'("# 2: protein.fasta")')
      write(4,'("# Matrix: EBLOSUM62")')
      write(4,'("# Gap_penalty: 10.0")')
      write(4,'("# Extend_penalty: 0.5")')
      write(4,'("#")')
      write(4,'("# Length: 432")')
      write(4,'("# Identity:     431/432 (99.8%)")')
      write(4,'("# Similarity:   431/432 (99.8%)")')
      write(4,'("# Gaps:           1/432 ( 0.2%)")')
      write(4,'("# Score: 2162.0")')
      write(4,'("#")')
      write(4,'("#")')
      write(4,'("#=======================================")')

      do i=1,nres_nascent
        seq3(i:i)='-'
      enddo

C NASCENT - blank off the residues after the final nascent one

      do i=nres_nascent+1,nres
        seq_of_fasta(i:i)=' '
      enddo

      if(iter.eq.1) then
        seq2(1:1)='|'
        seq3(1:1)=seq_of_fasta(1:1)
      elseif(iter.eq.2) then

C change to code 10/18/2021 - the existing code did not work with short
C nascent peptides - not sure why but there was a mismatch with the
C Cterminal residue being sometimes refered to as residue 'N' and other
C times as residue '1' - here I'm going to just make sure that all
C residues up to the final one are aligned - it may work...

C       seq2(nres_nascent:nres_nascent)='|'
C       seq3(nres_nascent:nres_nascent)=
C    &       seq_of_fasta(nres_nascent:nres_nascent)

        do j=1,nres_nascent
          seq2(j:j)='|'
        enddo
ccc     seq2(nres_nascent:nres_nascent)='|'
        seq3(1:nres_nascent)=seq_of_fasta(1:nres_nascent)

      endif

      iround=0
      do mp=1,nres_nascent,50
        iround=iround+1
        ibeg=(iround-1)*50+1
        iend=iround*50
        write(4,'("")')
        write(4,'("ENO_ECOLI      ",i5,1x,a50,i7)')ibeg,
     &             seq_of_fasta(ibeg:iend),iend
        write(4,'("               ",6x,a50)')seq2(ibeg:iend)
        write(4,'("protein.fasta  ",i5,1x,a50,i7)')ibeg,
     &             seq3(ibeg:iend),iend
      enddo

      write(4,'("")')
      write(4,'("")')
      write(4,'("#=======================================")')
      write(4,'("#=======================================")')

      stop
      end

         SUBROUTINE init_random_seed()

            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)

            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)

         END SUBROUTINE


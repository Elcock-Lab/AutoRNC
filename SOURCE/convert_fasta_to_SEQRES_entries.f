
C this is a preprocessing code that reads a .fasta file and generates a
C corresponding .seqres file

      program make_seqres

      integer iatom,ile,ilr,ires,ilast,ihave,num_trans
      integer num_mons_orig
      character*12 my_fmt
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
      character*70 outstring
      character*3  rnam
      character*80 filein,fileout,filetrans,crap
      real   x,y,z
      call getarg(1,crap)
      read(crap,*) filein  
      call getarg(2,crap)
      read(crap,*) fileout  
      call getarg(3,crap)
      read(crap,*) chainID  
199   format(a80)
c
      open(2,file=filein,status='old')
      open(3,file=fileout,status='unknown')

C read the fasta files for the input sequences

        nres=0
        read(2,*) ! skip top line
26      read(2,27,end=28)seq_chunk
27      format(a10000)
c27      format(a60)

        do i=1,len(trim(seq_chunk)) ! line don't need to have 60 chars
          if(seq_chunk(i:i).ne.' ') then
            nres=nres+1
            seq_of_fasta(nres:nres)=seq_chunk(i:i)
          else
            goto 28
          endif
        enddo

        goto 26
28      close(2)

        do i=1,70
          outstring(i:i)=' '
        enddo
        outstring(1:19)='SEQRES   1 X  123  '
        outstring(12:12)=chainID
        write(outstring(13:17),'(i5)')nres
        nline=0
        do mp=1,nres,13
          nline=nline+1
          write(outstring(7:10),'(i4)')nline
          ibeg=mp
          iend=min(nres,mp+12)
          j=0
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
          enddo

C remember to blank out the rest of the last line

          do j=i2+1,70
            outstring(j:j)=' '
          enddo

          write(3,'(a70)')outstring(1:70)
        enddo
      close(3)
      stop
      end

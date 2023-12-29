
C this is the main code used to generate a script for producing RNC
C models with AutoRNC - it has a *lot* of arguments, but once figured
C out it does the job pretty seamlessly
 
      program make_script_for_autornc
 
      implicit real(a-h,o-z)

      character*200 junk
      character*500 my_fmt2

      character*1   chainID
      character*3   char3
      character*4   char4
      character*6   char6
      character*8   char8
      character*9   char9
      character*10  char10
      character*20  char20A
      character*20  char20B
      character*24  char24
      character*26  char26
      character*27  char27
      character*29  char29
      character*52  char52
      character*58  char58
      character*62  char62
      character*72  char72
      character*80  char80
      character*83  char83
      character*137 char137
      character*137 char137B
      character*152 char152
      character*190 char190
      character*200 char200
      character*80  fname_TUs
      character*80  fname_footprints
      character*80  fname_gene2uniprot
      character*12  folder_name
      character*200 filebloks
      character*200 filesselems
      character*200 ribopath
      character*200 filepath
      character*200 filetmp
      character*200 codepath
      character*200 dipeppath
      character*32  filename

      integer       ibeg_blok(10000)
      integer       iend_blok(10000)
      parameter    (ngene_max=100)
      character*4   gene_nam(ngene_max)
      integer       gene_beg(ngene_max)
      integer       gene_end(ngene_max)
      integer       TU_beg
      integer       TU_end
      integer       TU_len

      parameter    (nribo_max=200)
      integer       ribo_beg(nribo_max)
      integer       ribo_end(nribo_max)
      integer       my_gene(nribo_max)
      character*4   my_gene_nam(nribo_max)
      character*6   my_uniprotID(nribo_max)
      integer       last_residue(nribo_max)

      character*6   uniprotID(10000)
      character*4   gene_of_uniprotID(10000)

      logical       rna_folded
      logical       protein_folded
      logical       iseparate
      logical       trROSETTA
      logical       AHEMODEL
      logical       oligomer
      logical       monomer
      logical       omit_TF

      logical       filebloks_exists
      logical       filesselems_exists
      logical       use_ribosome

      integer       sidechain_mode
      integer       nmodels

      call getarg(1,junk)
      read(junk,*)char6      ! jobID - uniprotID generally
      call getarg(2,junk)
      read(junk,*)nlength    ! length of construct in residues
      call getarg(3,junk)
      read(junk,*)nstart     ! first res to build in Nter direction
      call getarg(4,junk)
      read(junk,*)filebloks  ! file listing all folded blocks
      call getarg(5,junk)
      read(junk,*)filesselems! file listing all SS elements
      call getarg(6,junk)
      read(junk,*)ribopath   ! path to ribosome structure
      call getarg(7,junk)
      read(junk,*)filepath   ! path to protein structure 
      call getarg(8,junk)
      read(junk,*)codepath   ! path to AHE's code
      call getarg(9,junk)
      read(junk,*)dipeppath  ! 2023 path to dipeptide libraries
      call getarg(10,junk)
      read(junk,*)sidechain_mode ! =0 for all atoms, =2 for no SC atoms
      call getarg(11,junk)   ! 2023
      read(junk,*)nmodels    ! #models desired
      call getarg(12,junk)   ! 2023
      read(junk,*)char3      ! ="yes" to use ribosome
      call getarg(13,junk)   ! 2023
      read(junk,*)d_intra    ! intra-chain clash cutoff
      call getarg(14,junk)   ! 2023
      read(junk,*)d_inter    ! inter-chain clash cutoff
      call getarg(15,junk)   ! 2023
      read(junk,*)nback      ! nbacktrack
      call getarg(16,junk)   ! 2023
      read(junk,*)iblock     ! nascent chain res to block (if any)

C 2023 - set use_ribosome logical

      if(char3(1:3).eq.'yes'.or.char3(1:3).eq.'YES') then
        write(*,*)
        write(*,*)'arg11 = ',char3,' so we use ribosome'
        write(*,*)
        use_ribosome=.true.
      else
        write(*,*)
        write(*,*)'arg11 = ',char3,' so we do NOT use ribosome'
        write(*,*)'set it to yes or YES if you want it!'
        write(*,*)
        use_ribosome=.false.
      endif 

      call init_random_seed()

      AHEMODEL=.true.

      chainID="A" ! this should be our default value

C commenting out the references to oligomer here..

      oligomer=.false.
      monomer=.true. 

c     if(char10(1:8).ne.'oligomer'.and.
c    &   char10(1:7).ne.'monomer') then
c       write(*,*)
c       write(*,*)'final argument must be either oligomer_X or monomer'
c       write(*,*)'where X is the chain ID of the one to be nascent'
c       write(*,*)'quitting now :('
c       write(*,*)
c       stop
c     elseif(char10(1:8).eq.'oligomer') then
c       if(trROSETTA) then
c         write(*,*)
c         write(*,*)'can only select "oligomer" with "AHEMODEL"'
c         write(*,*)'since trROSETTA models are all monomeric'
c         write(*,*)'quitting now :('
c         write(*,*)
c         stop
c       endif
c       oligomer=.true.
c       chainID=char10(10:10)
c       write(*,*)
c       write(*,*)'we will attempt to build a nascent oligomer'
c       write(*,*)'we will use chain ',chainID,' as new nascent chain'
c       write(*,*)'from AHEMODEL.pdb file - all others will be made'
c       write(*,*)'into HETATM entries...'
c       write(*,*)
c     elseif(char10(1:7).eq.'monomer') then
c       monomer=.true.
c       write(*,*)
c       write(*,*)'we will attempt to build a nascent monomer - this'
c       write(*,*)'is the default situation'
c       write(*,*)
c     endif

C find the length of the filepath and codepath
C this is crazy - fortran won't read a '/' as a slash - it stops reading
C at that point - so we will need to read the path in like this:

C ?home?aelcock?CURRENT_CODES etc

C the code will then change all '?' symbols to '/' and will add in a '/'
C at the end of the path if there isn't one already provided...

C UPDATE - this turns out to not be necessary *if* the entire path is
C surrounded on both sides by backslash+single-quotes, i.e. \'
C the following code works in either case so keep it as is...

      do n=1,200
        if(ribopath(n:n).eq.' ') then
          kpath=n-1
          exit
        endif
      enddo
      do n=1,kpath
        if(ribopath(n:n).eq.'?') ribopath(n:n)='/'
      enddo

      if(ribopath(kpath:kpath).ne.'/') then
        kpath=kpath+1
        ribopath(kpath:kpath)='/'
      endif
      write(*,*)'check ribopath = "',ribopath(1:kpath),'"'

      do n=1,200
        if(filepath(n:n).eq.' ') then
          lpath=n-1
          exit
        endif
      enddo
      do n=1,lpath
        if(filepath(n:n).eq.'?') filepath(n:n)='/'
      enddo

      if(filepath(lpath:lpath).ne.'/') then
        lpath=lpath+1
        filepath(lpath:lpath)='/'
      endif
      write(*,*)'check filepath = "',filepath(1:lpath),'"'

      do n=1,200
        if(codepath(n:n).eq.' ') then
          mpath=n-1
          exit
        endif
      enddo
      do n=1,mpath
        if(codepath(n:n).eq.'?') codepath(n:n)='/'
      enddo

      if(codepath(mpath:mpath).ne.'/') then
        mpath=mpath+1
        codepath(mpath:mpath)='/'
      endif
      write(*,*)'check codepath = "',codepath(1:mpath),'"'

      do n=1,200
        if(dipeppath(n:n).eq.' ') then
          npath=n-1
          exit
        endif
      enddo
      do n=1,npath
        if(dipeppath(n:n).eq.'?') dipeppath(n:n)='/'
      enddo

      if(dipeppath(npath:npath).ne.'/') then
        npath=npath+1
        dipeppath(npath:npath)='/'
      endif
      write(*,*)'check dipeppath = "',dipeppath(1:npath),'"'
      write(*,*)

C set some logicals and initial values here...

      protein_folded=.true. 
      nribos=1
      n=nribos
      my_uniprotID(nribos)=char6
 
C also set a logical if filesselems exists...

      filesselems_exists=.false.
      inquire(file=filesselems,exist=filesselems_exists)

C open up the comfile that we want to write all commands to

      open(unit=11,file='comfile_build_nascent_chain',
     &     status='unknown')

C first make the appropriate directory

      write(11,*) 
      write(11,'("# make directory and change to it")')
      write(11,*) 
      folder_name='AHEMODEL_XXX'
      write(char3,'(i3)')n
      if(char3(1:1).eq.' ') char3(1:1)='0'
      if(char3(2:2).eq.' ') char3(2:2)='0'
      if(char3(3:3).eq.' ') char3(3:3)='0'
      folder_name(10:12)=char3
      write(11,'("mkdir ",a12)')folder_name

C change to that directory

      write(11,'("cd ",a12)')folder_name

C now copy necessary files using the filepath provided on command line

      write(11,*) 
      write(11,'("# copy needed files from user-provided filepath")')
      write(11,*) 
      write(11,'("# RIBOSOME.pdb")')
      write(11,'("# NASCENT_CTERMINUS.pdb")')
      write(11,'("# input_file_for_ahemodel_2023_v1.0")')
      write(11,*) 

      write(my_fmt2,'(a,i0,a)'),'("cp ",a',kpath,
     & '"RIBOSOME.pdb .")'
      write(11,my_fmt2)ribopath(1:kpath)
      write(11,*) 

      write(my_fmt2,'(a,i0,a)'),'("cp ",a',kpath,
     & '"NASCENT_CTERMINUS.pdb .")'
      write(11,my_fmt2)ribopath(1:kpath)
      write(11,*) 

      write(my_fmt2,'(a,i0,a)'),'("cp ",a',kpath,
     &                           '"input_file_for_ahemodel_2023_v1.0_",
     &                    "NASCENT input_file_for_ahemodel_2023_v1.0")'
      write(11,my_fmt2)ribopath(1:kpath)

C if sidechain_mode<>0 then we need to remove sidechains from the
C RIBOSOME.pdb and NASCENT_CTERMINUS.pdb before proceeding...

      if(sidechain_mode.eq.1.or.sidechain_mode.eq.2) then

        write(11,*) 
        write(11,'("# remove appropriate sidechain atoms")')
        write(11,*) 

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &     '"remove_sidechains_in_pdb.exe RIBOSOME.pdb tmp.pdb ",i6)'
        write(11,my_fmt2)codepath(1:mpath),sidechain_mode
        write(11,'("mv tmp.pdb RIBOSOME.pdb")')
        write(11,*) 

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &     '"remove_sidechains_in_pdb.exe NASCENT_CTERMINUS.pdb ",
     &      "tmp.pdb ",i6)'
        write(11,my_fmt2)codepath(1:mpath),sidechain_mode
        write(11,'("mv tmp.pdb NASCENT_CTERMINUS.pdb")')

      endif

C now that we have manipulated RIBOSOME.pdb and NASCENT_CTERMINUS.pdb we
C can, if iblock>0, add blocking atoms around the specified residue...

      if(iblock.gt.0) then

        write(11,*) 
        write(11,'("# add in blocking atoms around nascent chain")')
        write(11,*) 

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &     '"add_blocking_atoms_to_ribosome_pdb.exe ",
     &      "NASCENT_CTERMINUS.pdb RIBOSOME.pdb BLOCKERS.pdb ",
     &      "20.0 3.0 ",i4)'
        write(11,my_fmt2)codepath(1:mpath),iblock
        write(11,'("cat BLOCKERS.pdb >> RIBOSOME.pdb")')

      endif

C execute a sed command to change the uniprotID - note that the
C substitution bit does NOT need to be in quotes (don't know why I always
C thought that it did...)

      write(11,*) 
      write(11,'("# set the uniprotID in the ahemodel input file")')
      write(11,'("# note that this must have rebuild_tails = yes")')
      write(11,*) 
      write(11,'("sed -i s/PXXXXX/",a6,
     &           "/ input_file_for_ahemodel_2023_v1.0")')my_uniprotID(n)

C execute a sed command to put in the sidechain_mode integer...

      write(11,*) 
      write(11,'("# set sidechain_mode in the ahemodel input file")')
      write(11,'("# =0 for all atoms; =1 for CBeta only; =2 for none")')
      write(11,*) 
      write(11,'("sed -i s/''SCMODE''/''",i6,
     &      "''/ input_file_for_ahemodel_2023_v1.0")')sidechain_mode

C 2023 execute a sed command to put in the nmodels...

      write(11,*)
      write(11,'("# set nmodels in the ahemodel input file")')
      write(11,*)
      write(11,'("sed -i s/''NMODELS''/''",i6,
     &      "''/ input_file_for_ahemodel_2023_v1.0")')nmodels

C 2023 execute a sed command to put in the intra-chain clash distance

      write(11,*)
      write(11,'("# set d_intra in the ahemodel input file")')
      write(11,*)
      write(11,'("sed -i s/''D_INTRA''/''",f6.2,
     &      "''/ input_file_for_ahemodel_2023_v1.0")')d_intra

C 2023 execute a sed command to put in the inter-chain clash distance

      write(11,*)
      write(11,'("# set d_inter in the ahemodel input file")')
      write(11,*)
      write(11,'("sed -i s/''D_INTER''/''",f6.2,
     &      "''/ input_file_for_ahemodel_2023_v1.0")')d_inter

C 2023 execute a sed command to put in nbacktrack

      write(11,*)
      write(11,'("# set nback   in the ahemodel input file")')
      write(11,*)
      write(11,'("sed -i s/''NBACK''/''",i5,
     &      "''/ input_file_for_ahemodel_2023_v1.0")')nback

C 2023 execute a sed command to put in codepath into the AHEMODEL input
C file - this should allow us to have all codes for AutoRNC and AHEMODEL
C present in the same folder...
C note that we use the separator '#' so that sed interprets slashes
C correctly for the path...

          write(my_fmt2,'(a,i0,a)'),
     &                '("sed -i s#''CODE_PATH''#''",a',mpath,'
     &      "''# input_file_for_ahemodel_2023_v1.0")'
          write(11,my_fmt2)codepath(1:mpath)

C 2023 same thing for dipeppath...

          write(my_fmt2,'(a,i0,a)'),
     &                '("sed -i s#''DIPEP_PATH''#''",a',npath,'
     &      "''# input_file_for_ahemodel_2023_v1.0")'
          write(11,my_fmt2)dipeppath(1:npath)

C now get the full-length .fasta

        write(11,*) 
        write(11,'("# get the construct .fasta from filepath")')
        write(11,*) 
        write(my_fmt2,'(a,i0,a)'),'("cp ",a',lpath,
     &                          'a6,".fasta .")'
        write(11,my_fmt2)filepath(1:lpath),my_uniprotID(n)

C now if filesselems_exists we want to make our fake .ss2 file (using
C the .fasta that we've just copied over) and we want to edit the input
C file for ahemodel to make sure it uses it...
C TEMP TEMP TEMP THIS NOT FINISHED YET

        if(filesselems_exists) then

          write(11,*) 
          write(11,'("# copy the list_of_ss_elems file and")')
          write(11,'("# use it to make a .ss2 file")')
          write(11,*) 

C pretty crude code here - we want a version of filesselems that has a
C generic name as it makes the call to the code much easier - so we just
C do a copy command here to generate PXXXXX.sselems

          do l=1,200
            if(filesselems(l:l).eq.' ') then
              ltmp=l-1
              exit
            endif
          enddo

C note that since filepath and filesselems are both of variable length
C it's easier to just concatenate them together and use that to make the
C writing of the copy statement doable - thanks Fortran!

          filetmp(1:lpath)=filepath(1:lpath)
          filetmp(lpath+1:lpath+ltmp)=filesselems(1:ltmp)

          write(my_fmt2,'(a,i0,a)'),'("cp ",a',lpath+ltmp,',
     &                                " ",a6,".sselems")'
          write(11,my_fmt2)filetmp(1:lpath+ltmp),my_uniprotID(n)

C now use that copy to make a .ss2 file to be read by AHEMODEL

          write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"make_fake_ss2_file.exe ",a6,".sselems ",a6,".fasta ",
     &         i6," ",a6,".ss2")'
          write(11,my_fmt2)codepath(1:mpath),my_uniprotID(n),
     &                     my_uniprotID(n),nlength,my_uniprotID(n)

C finally, make sure that the AHEMODEL input file looks for it...
C note that this assumes that the starting input file says "junk.ss2"
C also note that use of the '' to put in a single ' - this is the same
C trick as used for double quotes...

          write(11,'("sed -i s/''junk.ss2  ''/''",a6,
     &   ".ss2''/ input_file_for_ahemodel_2023_v1.0")')my_uniprotID(n)
                  
        endif

C now make an AHEMODEL template for the C-terminal most residue

        write(11,*) 
        write(11,'("# start template construction with C-ter chunk")')
        write(11,*) 

C arguments:

C [1] input .fasta
C [2] input pdb file for C-terminal most residue of nascent chain
C [3] output sequence alignment file for AHEMODEL
C [4] output TEMPLATE pdb file for AHEMODEL
C [5] input C-terminal-most residue of construct (i.e. length)
C [6] first residue to build in N-ter direction

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"prepare_ahemodel_for_folded_nascent_protein_final.exe ",
     &         a6,".fasta NASCENT_CTERMINUS.pdb ",a6,"_TEMPLATE_",
     &         "ahemodel.001.out TEMPLATE_ahemodel.CTERMINUS.pdb ",
     &         2i6)'
        write(11,my_fmt2)codepath(1:mpath),my_uniprotID(n),
     &                   my_uniprotID(n),nlength,nstart

C now get the full-length AlphaFold2.pdb or AHEMODEL.pdb etc

        write(11,*) 
        write(11,'("# continue template using structural model")')
        write(11,*) 
        write(my_fmt2,'(a,i0,a)'),'("cp ",a',lpath,
     &                             'a6,".pdb model0.pdb")'
        write(11,my_fmt2)filepath(1:lpath),my_uniprotID(n)

C get all atoms that belong to the desired chain here
C the format here is such a pain in the ass to figure out that we'll
C just dump all the text into a string and then write that out

        char200(  1: 30)='grep "^.\{21\}?" model0.pdb | '
        char200( 31: 55)='grep "^ATOM" > model1.pdb'
        char200( 15: 15)=chainID
        write(11,'(a55)')char200(1:55)

C now we want to remove any hydrogens...

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"remove_hydrogens_in_pdb.exe model1.pdb model2.pdb")'
        write(11,my_fmt2)codepath(1:mpath)

C finally, we may want to remove sidechain atoms...

        if(sidechain_mode.eq.0) then

          write(11,'("cp model2.pdb model3.pdb")')

        elseif(sidechain_mode.eq.1) then

          write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &      '"remove_sidechains_in_pdb.exe model2.pdb model3.pdb 1")'
          write(11,my_fmt2)codepath(1:mpath)

        elseif(sidechain_mode.eq.2) then

          write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &      '"remove_sidechains_in_pdb.exe model2.pdb model3.pdb 2")'
          write(11,my_fmt2)codepath(1:mpath)

        endif

C and use it to make a full-length template, from which folded chunks
C can be extracted as desired...

        write(11,*) 
        write(11,'("# make template from provided structure")')
        write(11,*) 

        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"convert_fasta_to_SEQRES_entries.exe ",a6,
     &        ".fasta model3.seqres A")'
        write(11,my_fmt2)codepath(1:mpath),my_uniprotID(n)

        write(11,'("cat model3.seqres model3.pdb > model3_final.pdb")')

C if we're using an AHEMODEL structure we need to modify the
C occupancy fields on any bits that were modeled in...

        write(11,*) 
        write(11,'("# if AHEMOD struct then modify occupancy fields")')
        write(11,*) 
        write(11,'("sed -i s/""  0.00  0.00""/""  1.00  0.00""/ ",
     &             "model3_final.pdb")')
        write(11,'("sed -i s/""  1.00 99.99""/""  1.00  0.00""/ ",
     &             "model3_final.pdb")')

        write(11,*) 
        write(11,'("# do a final clean-up of the template")')
        write(11,*) 
        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"clean_up_template.exe model3_final.pdb ",
     &        "TEMPLATE_ahemodel.FULL_LENGTH.pdb XXXX")'
        write(11,my_fmt2)codepath(1:mpath)

C now we need to open up the filebloks file and find out how many (if
C any folded blocks) we need to extract

        nbloks=0
        filebloks_exists=.false.
        inquire(file=filebloks,exist=filebloks_exists)

        if(filebloks_exists) then
          open(unit=15,file=filebloks,status='unknown')
          read(15,'(1x)') ! skip title line
801       read(15,'(a80)',end=802)char80
          read(char80,*)i,j,k
          nbloks=nbloks+1
          ibeg_blok(nbloks)=j
          iend_blok(nbloks)=k
          goto 801
802       close(15)
        endif

C now cut out only those structured residues from the native state pdb

        write(11,*) 
        write(11,'("# cut out the structured residues from the")')
        write(11,'("# full-length template into a culled version")')
        write(11,*) 
        do nb=1,nbloks
          write(char6,'(i6)')nb
          do m=1,5
            if(char6(m:m).eq.' ') char6(m:m)='0'
          enddo
          write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"cull_full_length_template_by_folding_blok.exe ",
     &        "TEMPLATE_ahemodel.FULL_LENGTH.pdb ",
     &        "TEMPLATE_ahemodel.CULLED.pdb.",a6,3i6)'
          write(11,my_fmt2)codepath(1:mpath),char6,ibeg_blok(nb),
     &                                    iend_blok(nb),nlength
        enddo

C now if "oligomer" was selected *and* it's an AHEMODEL case then add on
C the non-chain A atoms as HETATM entries here...

c       if(oligomer) then
c         write(11,*) 
c         write(11,'("# since we selected oligomer we need to get")')
c         write(11,'("# the other chains and make them HETATMs")')
c         write(11,*) 
c         char200(  1: 33)='grep -v "^.\{21\}?" model0.pdb | '
c         char200( 34: 55)='grep "^ATOM" > tmp.pdb'
c         char200( 18: 18)=chainID
c         write(11,'(a55)')char200(1:55)
c         write(11,'("sed -i s/""ATOM  ""/""HETATM""/ tmp.pdb")')
c         write(11,'("cat tmp.pdb >> TEMPLATE_ahemodel.CULLED.pdb")')
c       endif

C now combine the templates, after listing them in list_of_pdbs
C note the neat trick being used here: when "" appears it is interpreted
C as "put in a single double-quote" - i.e. it puts in a single " ...
C this saves a *lot* of messing around with single and double quotes...
C that idea was taken from here:

C http://computer-programming-forum.com/49-fortran/213d85357b013f7f.htm

        write(11,*) 
        write(11,'("# prepare to combine our truncated nascent chain")')
        write(11,'("# template with our Cterm residue")')
        write(11,*) 

        write(11,'("rm list_of_pdbs")')
        do nb=1,nbloks
          write(char6,'(i6)')nb
          do m=1,5
            if(char6(m:m).eq.' ') char6(m:m)='0'
          enddo
          write(11,'("echo ""TEMPLATE_ahemodel.CULLED.pdb.",a6,
     &               "    1"" >> list_of_pdbs")')char6
        enddo
        write(11,'("echo ""TEMPLATE_ahemodel.CTERMINUS.pdb",
     &             "        1"" >> list_of_pdbs")')

C arguments:

C [1] list of templates
C [2] output TEMPLATE_ahemodel.COMBINED.pdb
C [3] =1 since we just want to use SEQRES from the first pdb file (we do:
C        it's the only one that's complete)
C [4] =0 since we don't want to clip off the first and last residues of
C        each structured domain
C [5] =2 since we want the 2nd pdb file to be the one that doesn't move

        write(11,*) 
        write(11,'("# combine the two with the following command")')
        write(11,'("# 1 means take SEQRES from first listed pdb")')
        write(11,'("# 0 means do not clip residues from domains")')
        write(11,'("# I means use Ith pdb file as stationary - this")')
        write(11,'("# will be the C-terminal chunk as this is")')
        write(11,'("# placed in the ribosome structure...")')
        write(11,*) 
        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &       '"combine_templates_handles_one_oligomer.exe ",
     &        "list_of_pdbs TEMPLATE_ahemodel.COMBINED.pdb 1 0 ",i6)'
        write(11,my_fmt2)codepath(1:mpath),nbloks+1

C now add on the ribosome atoms:
C note that we have two ways of doing this depending on whether we're
C choosing to omit TF from the ribosome...
C 2023 only do this if use_ribosome

        if(use_ribosome) then
          write(11,*) 
          write(11,'("# cat on the ribosome HETATMs to the template")')
          write(11,*) 
          write(11,'("cat RIBOSOME.pdb >> ",
     &                 "TEMPLATE_ahemodel.COMBINED.pdb")')
        endif

C 2021 for each of the following lines we remove advance="no" *AND* we
C remove the semicolon at the end of the line...
 
        write(11,*) 
        write(11,'("# copy the final combined template pdb file")')
        write(11,'("# to the generic name that AHEMODEL looks for")')
        write(11,*) 
        write(11,'("cp TEMPLATE_ahemodel.COMBINED.pdb ",
     &             "TEMPLATE_ahemodel.pdb ")')

C note that we need to be told the first residue to build moving in the
C N-ter direction on the command line - this will be the first residue
C prior to the first residue of the ribosome-bound C-terminus...

        write(11,*) 
        write(11,'("# make tail_rebuilds.txt")')
        write(11,*) 
        write(11,'("echo "" N   1 ",i5,""" > tail_rebuilds.txt ")')
     &             nstart

        write(11,*) 
        write(11,'("# run AHEMODEL!!!")')
        write(11,*) 
        write(my_fmt2,'(a,i0,a)'),'(a',mpath,
     &  '"ahemodel_2023_v1.2.exe < input_file_for_ahemodel_2023_v1.0")'
        write(11,my_fmt2)codepath(1:mpath)

C make a FINAL copy of the protein before going on

C       write(11,*) 
C       write(11,'("# make a FINAL copy of the AHEMODEL output")')
C       write(11,*) 
C       write(11,'("cp ",a6,".AHEMODEL.pdb ",
C    &             a6,".AHEMODEL.FINAL.pdb ")')
C    &             my_uniprotID(n),my_uniprotID(n)

C close the comfile

        close(11)

      call system("chmod +x comfile*")

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


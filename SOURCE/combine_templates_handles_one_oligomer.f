
C this is preprocessing code that is used to combine a series of
C template fragments into a single template that is provided to
C ahemodel_2023_v1.2.exe, which is the basic engine for AutoRNC
C note: users of AutoRNC should not need to intervene in this code

C code UPDATED to allow very long linkers (>1000 residues) - used to
C allow only 100 residues (changed char100 and char90 to char1000 and
C char990 respectively - also set nlinker_seq to have 990 characters)

C this code reads template.pdb files and merges them

      program place_molecules
c
      implicit real(a-h,o-z)
      parameter   (pi=3.141592654)
      INTEGER, PARAMETER :: K4B=selected_int_kind(9)

      interface
      subroutine pdbsup(nat_loc,rm_loc,rf_loc)

      integer, intent(in)     :: nat_loc
      real, intent(inout)     :: rm_loc(1:nat_loc,1:3)
      real, intent(in)        :: rf_loc(1:3,1:3)

      end subroutine
      end interface


      integer      :: omp_get_thread_num

      character*1 atail
      character*4 atmnam   
      character*3 resnam 
      character*89 string_sup
      integer     resnum,num_box

      real, allocatable :: af(:)
      real, allocatable :: bf(:)
      real, allocatable :: cf(:)
c     real, allocatable :: xf(:)
c     real, allocatable :: yf(:)
c     real, allocatable :: zf(:)
      real, allocatable :: tx(:)
      real, allocatable :: ty(:)
      real, allocatable :: tz(:)
      real, allocatable :: ax1(:)
      real, allocatable :: ax2(:)
      real, allocatable :: ax3(:)
      real, allocatable :: ay1(:)
      real, allocatable :: ay2(:)
      real, allocatable :: ay3(:)
      real, allocatable :: az1(:)
      real, allocatable :: az2(:)
      real, allocatable :: az3(:)

      real, allocatable :: xacc(:,:)
      real, allocatable :: yacc(:,:)
      real, allocatable :: zacc(:,:)
      real, allocatable :: rf(:,:)
      real, allocatable :: rm(:,:)
      real                :: cm_pdbsup(3)
      real                :: cf_pdbsup(3)
      real                :: t_pdbsup(3,3)

      type atomsplaced
        real                                  :: i
        real                                  :: j
        real                                  :: k
        real                                  :: x
        real                                  :: y
        real                                  :: z
        integer                               :: n ! residue number
        integer                               :: m ! molecule number
        integer                               :: t ! molecule type   
        character*1                           :: c ! chain identifier
        character*4                           :: a ! atom name         
        character*3                           :: r ! residue name       
      end type atomsplaced
      type (atomsplaced), allocatable :: coo_arr(:)

      type moltypes
        real        ,allocatable              :: xs(:) ! x of src
        real        ,allocatable              :: ys(:) ! y of src
        real        ,allocatable              :: zs(:) ! z of src
        integer     ,allocatable              :: ts(:) ! typ of src
        real        ,allocatable              :: rs(:) ! rate decay
        real        ,allocatable              :: x(:)
        real        ,allocatable              :: y(:)
        real        ,allocatable              :: z(:)
        integer     ,allocatable              :: n(:) ! residue number
        character*1 ,allocatable              :: c(:) ! chain identifier
        character*4 ,allocatable              :: a(:) ! atom name         
        character*3 ,allocatable              :: r(:) ! residue name       
        integer                               :: num_req
        integer                               :: natoms
        integer                               :: add_typ
        integer                               :: add_src
        integer                               :: add_sum
      end type moltypes
      type (moltypes), allocatable  :: addeds(:)

      type gridtype
        integer                     :: num
        real, dimension(:), pointer :: x
        real, dimension(:), pointer :: y
        real, dimension(:), pointer :: z
      end type gridtype
      type (gridtype), allocatable  :: grid (:,:,:)                                                                                
      type (gridtype)  :: temp_arr                                                                                                 

      real, allocatable :: prob_grid(:,:,:)
     
      integer,     allocatable :: nname_seq(:)
      character*80,allocatable :: fname_seq(:)
      character*80,allocatable :: fname_pdb(:)
      character*80,allocatable :: fname_pdb_chn(:,:)
      character*80,allocatable :: fasta_pdb_chn(:,:)
      integer     ,allocatable :: num_chns(:)
      character*10000,allocatable :: seq_of_fasta(:)
      character*60                :: seq_chunk
      character*10000,allocatable :: seq_of_pdb_chn(:,:)
      integer,        allocatable :: iam_in_pdb_chn(:,:,:)
      integer,        allocatable :: iam_in_SS_elem(:,:,:)
      character*10000             :: seq_scwrl
      integer        ,allocatable :: len_of_pdb_chn(:,:)
      integer        ,allocatable :: num_res_in_seq(:)
      character*1                 :: chn_last
      character*22                :: loopy_file

      integer num_res_in_frag(100000)
      integer idone(10000)
      integer ibeg_word(10000)
      integer iend_word(10000)
      integer ibeg_loop(10000)
      integer iend_loop(10000)
      integer ibeg_loop_final(10000)
      integer iend_loop_final(10000)
      integer i_not_made_yet(10000)
      integer i_in_pdb_fasta(10000) ! only for res in current chain
      integer i_in_pdb_fasta_tot(10000) ! for all res in all chains
      integer i_in_pdb_coord(10000)
      integer i_am_in_a_del(10000)
      integer my_nnn(10000)
      integer my_mmm(10000)
      integer my_nres2(10000)

      integer i_old_unstruc(100000) ! =1 if unstruc in original pdb
      integer i_am_original(100000) ! =1 if present in original pdb
      integer i_am_by_loopy(100000) ! =1 if made by loopy
      integer num_tot_res
      integer num_tot_nters
      integer num_tot_cters
      integer num_tot_loops
      integer my_frst_res(100000) ! dimensioned 1:# total chains
      integer my_last_res(100000)
      integer my_totl_res(100000)
      integer irebuild_nters(100000)
      integer irebuild_cters(100000)
      integer ifrg_tot_nters(100000)
      integer ifrg_tot_cters(100000)
      integer ifrg_tot_loops(100000)
      integer ibeg_tot_nters(100000)
      integer iend_tot_nters(100000)
      integer ibeg_tot_cters(100000)
      integer iend_tot_cters(100000)
      integer ibeg_tot_loops(100000)
      integer iend_tot_loops(100000)
      integer ibeg_loc_nters(100000)
      integer iend_loc_nters(100000)
      integer ibeg_loc_cters(100000)
      integer iend_loc_cters(100000)
      integer ibeg_loc_loops(100000)
      integer iend_loc_loops(100000)
      integer num_act_cters          ! # of actual already-placed Cters
      integer ipos_act_cters(100000) ! scwrl res# of the above
      character*1 ichn_tot_nters(100000)
      character*1 ichn_tot_cters(100000)
      character*1 ichn_tot_loops(100000)
      character*1000 char1000
      character*990 char990
      character*990 nlinker_seq(1000)
      character*66 char66
      character*6  char6
      character*8  char8
      integer my_scwrl_res(100000)
      integer my_tot_res(100000)
      integer i_done_res(100000)
      real rand_res(1000000)

      integer match_m1(1000,1000)
      integer match_m2(1000,1000)
      integer ioffset(1000)

      integer my_chains_seqres_entry(1000)

      integer my_seqres_chain(1000)
      integer my_seqres_chain_type(1000)
      integer my_seqres_chain_type_done(1000)
      integer              :: nres_seqres(1000)   
      character*3          :: rnam_seqres(1000,100000)   
      integer              :: nres_seqres_type(1000)   
      character*3          :: rnam_seqres_type(1000,100000)   

      integer my_atom_chain(1000)
      integer my_atom_chain_type(1000)
      integer my_atom_chain_type_done(1000)
      integer ibeg_seg(1000)
      integer iend_seg(1000)
      integer indx(1000)
      character*3 rnxm(1000)
      character*1 chain_2find

      real a(3),b(3),c(3),vca(3),vcb(3),xp(3),yp(3),zp(3)

      logical use_SS_constraints
      logical rebuild_tails

      real xntm(1000000)
      real yntm(1000000)
      real zntm(1000000)
      real xctm(1000000)
      real yctm(1000000)
      real zctm(1000000)
      integer jf(1000000) ! res# of amino acids already there
      integer jl(1000000) ! res# of amino acids used for clash checking
      real xf(1000000) ! amino acids already there
      real yf(1000000)
      real zf(1000000)
      real xh(1000000) ! HETATM entries already there
      real yh(1000000)
      real zh(1000000)
      real xm(10000) ! amino acid to be added
      real ym(10000)
      real zm(10000)
c     real xc(10000) ! amino acid adjacent to one to be added
c     real yc(10000)
c     real zc(10000)
      real xl(10000) ! amino acids close to the one to be added
      real yl(10000)
      real zl(10000)
      real xll(10000) ! used for checking whether loopy loops are closed
      real yll(10000)
      real zll(10000)
      character*100000 seq_tot_res
      character*100000 chn_tot_res
      integer          res_tot_res(100000)
      character*12 my_fmt
      character*80 my_fmt2
      character*3  ajunk,res_n
      character*1  chain
      character*80 final_name
      character*80 fname_inp 
      character*80 fname_mol 
      character*80 fname_out 
      integer      fname_len
      character*80 atemp       
      character*80 junk        
      character*4  addtype     
      character*80 string
      character*80 strang(10000) ! ATOMS
      character*80 streng(10000) ! HETATMS
      character*15    temp4_nter
      character*15    temp4_cter
      character*35    sys_stringy
      character*77    sys_stringx
      character*1000  sys_string1
      character*1000  sys_string2
      character*1000  sys_string3
      character*1000  sys_string4
      character*1000  sys_string5
      character*1000  sys_string6
      character*1000  sys_string7
      character*1000  sys_string8
      character*1000  sys_string9
      character*80    new_string(1:1000)
      character*80    strong
      character*80    last_strong
      character*100 strung1
      character*100 strung2
      character*100 strung3
      character*100 strung4
      character*100 strung5
      character*90 streng0
      character*90 streng1
      character*90 streng2
      character*4  num_of_res
      character*1  cc
      character*50 a1,a2,a3,a4,a5,a6,a7,a8,a9
      character*50,allocatable :: falign(:,:,:)
      character*50,allocatable :: falign2(:,:,:) ! for needle rewrite
      integer,     allocatable :: ialign(:,:,:)
      real,        allocatable :: salign(:,:,:)
      character*78 blah
      character*1  char1
      character*20 dipep
      character*2  dipep_cur
      character*3  nam_het_typs(100000)
      integer num_peps(1:20,1:20)
      integer align_mode
      integer membrane_dive

      character*1 last_ss_type
      character*1  curr_ss_type
      character*1  my_dssp_type(8)
      integer      nss(1000000)
      integer      do_me_still(1000000)
      integer      ibeg_ss(1000000)
      integer      iend_ss(0:1000000)
      integer      ilen_ss(1000000) ! length of SS element
      integer      ibeg_ss_tmp(1000000)
      integer      iend_ss_tmp(1000000)
 
      integer      tmplte_restyp(1000000)
      integer      target_restyp(1000000)
    
      character*1          :: sequenceA(100000)
      character*1          :: sequenceB(100000)
      character*1          :: sequenceM(100000)
      character*1          :: sequenceS(100000) ! SS element
      character*1          :: sequenceT(100000) ! in pdb?
      character*1          :: sequenceU(100000) ! '0123...'
      integer, allocatable :: seq1(:) 
      integer, allocatable :: seq2(:) 
      real,    allocatable :: score(:,:) 
      real,    allocatable :: val(:,:) 
      integer, allocatable :: idir(:,:) 
      real,    allocatable :: jpV(:,:) 
      real,    allocatable :: jpH(:,:) 
      real,    allocatable :: preV(:,:) 
      real,    allocatable :: preH(:,:) 
      real,    allocatable :: j2i(:)

C allowed dimensions here are complete overkill:

C up to 1000 chains
C up to 10000 residues in each chain

      integer              :: iwantthis(10000) ! up to 1000 atms in res
      character*4          :: atom_name(10000) ! up to 1000 atms in res
      real                 :: occupancy(10000) ! up to 1000 atms in res
      integer              :: ires_real(10000) ! up to 1000 atms in res
      real                 :: occ_A,occ_B
      character*70         :: outstring
      integer              :: rnum,rnum_frst,rnum_last,rnum_orig
      character*4          :: anam
      character*3          :: rnam
      character*1          :: cnam,cnam_last,cnam_output
      character*1          :: dnam,dnam_last,dnam_output
      character*1          :: fasta
      character*10000      :: seq_seqres(1000)   
      character*10000      :: seq_pisa(1000)   
      integer              :: my_real_rnum_pisa(1000,10000)   
      integer              :: nres_pisa(1000)   
      character*100        :: seqres_pdb
      character*100        :: initial_pdb(1000)
      character*100        :: final_pdb
      character*100        :: filename
      character*100        :: list_of_pdbs
      character*1          :: chain_seqres(1000)
      character*1          :: chain_pisa(1000)
      character*1          :: mychain_pdb(1000)
      integer              :: IDchain_pdb(1000)
      integer              :: nchains_pdb(1000)
      integer              :: nlinker_res(1000)
      integer              :: nres_first(1000)
      integer              :: my_match(1000)
      integer              :: my_offset(1000)
      integer              :: rnum_offset
      character*1          :: res_pisa
      character*1          :: res_seqres

      real                 :: xcen(1000)
      real                 :: ycen(1000)
      real                 :: zcen(1000)
      character*66         :: string66
      character*66         :: string66_stored(10000)
      character*4          :: char4,pdbid
      character*7          :: char7
      character*1          :: chain_old
      character*1          :: chain_new
      character*1          :: chainID(10000)
      character*1          :: chain_stored(10000)
      character*13         :: temp_pdb
      character*14         :: temp_pdb2
      character*10         :: char10
      character*10         :: char10_stored(10000)
      character*6          :: occpy
      character*3          :: res_nam_final(100000)

      integer              :: chain_is_too_short(10000)
      integer              :: chain_is_nucleic(10000)
      real                 :: xnf(100000),ynf(100000),znf(100000)
      real                 :: xcf(100000),ycf(100000),zcf(100000)
      integer              :: my_res_final(100000)
      integer              :: iskip(100000)
      integer              :: frst_line(100000)
      integer              :: last_line(100000)
      integer              :: ifoundN(100000)
      integer              :: ifoundCA(100000)
      integer              :: ifoundC(100000)
      integer              :: ifoundO(100000)

      integer              :: skipfirstlast

      character*1000       :: strong1000 
      character*1000       :: string1000 
      character*3          :: acrap1,acrap2,acrap3,acrap4
      character*3          :: mod_res(511)
      character*3          :: ori_res(511)
      character*1 aa(20)
      data aa /
     &   'A','R','N','D','C','Q','E','G','H','I',
     &   'L','K','M','F','P','S','T','W','Y','V' /

      blah(1:26)='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      blah(27:52)='abcdefghijklmnopqrstuvwxyz'
      blah(53:78)='!@#$%^&*()_+-={}[]|\:;<>,.'

      dipep(1:20)='ACDEFGHIKLMNPQRSTVWY'

      write(*,*)
      write(*,*)'if the code fails on reading the arguments'
      write(*,*)'you may need to add a 3rd arg =1'
      write(*,*)'and a 4th arg=1 generally'
      write(*,*)

      call getarg(1,junk)
      read(junk,*)list_of_pdbs ! original pdb file with alt confs
      call getarg(2,junk)
      read(junk,*)final_pdb  ! output pdb file witout alt confs
      call getarg(3,junk)
      read(junk,*)ialready   ! =1 if all pdb files already have the 
                             ! same SEQRES entries - no need to add
                             ! linkers and no need to read any SEQRES
                             ! entries for any pdbs other than 1
      call getarg(4,junk)
      read(junk,*)skipfirstlast ! =1 if we skip the first and last
                                ! residues of each domain that links
                                ! to another
      call getarg(5,junk)
      read(junk,*)nonmovingpdb  ! if<1 then determine automatically

C each line of list_of_pdbs contains pdb name, then # of residues that
C make up the linker, then the sequence to be assigned to it - the
C latter can be X if necessary...

      open(unit=11,file=list_of_pdbs,status='unknown')
      open(unit=12,file=final_pdb,status='unknown')

C initial_pdb(:) is the list of pdb files
C nlinker_res(:) is the number of linker res needed *after* this pdb
C nlinker_seq(:) is the one-letter-code sequence of the linker

      num_pdbs=0
10    read(11,'(a1000)',end=20)char1000
      read(char1000,*)filename,nlink
      num_pdbs=num_pdbs+1
      initial_pdb(num_pdbs)=filename
      if(ialready.eq.1) then 
        write(*,*)
        write(*,*)'linker entries are ignored since ialready=1 '
        write(*,*)
      else  
        nlinker_res(num_pdbs)=nlink
        if(nlink.gt.0) read(char1000,*)filename,nlink2,char990(1:nlink)
        nlinker_seq(num_pdbs)(1:nlinker_res(num_pdbs))=char990(1:nlink)
      endif

      goto 10
20    close(11)
 
      write(*,*)'number of template pdbs = ',num_pdbs

C now open all of them and write all non-ATOM/non-HETA entries and
C non-SEQRES entries to the final_pdb - these *should* just be REMARK
C entries added to the templates so we know their provenance...

      do n=1,num_pdbs
        open(unit=20+n,file=initial_pdb(n),status='unknown')
30      read(20+n,'(a66)',end=40)char66
        if(char66(1:4).eq.'ATOM') goto 30
        if(char66(1:4).eq.'HETA') goto 30
        if(char66(1:4).eq.'SEQR') goto 30
        write(12,'(a66)')char66
        goto 30
40      close(20+n)
      enddo

C now, before doing anything else, we'll read the initial_pdb files,
C look at the SEQRES entries and determine how many chains are in each
C at this stage of code-development we allow only one of the initial_pdb
C files to have more than one chain...
C we will also store the chainID of the first chain so we know when to
C stop reading the actual SEQRES entries later...

      nchains_pdb=0
      nchain_max=0
      noligomers=0

      do n=1,num_pdbs

        open(unit=20+n,file=initial_pdb(n),status='unknown')

        nchain_seqres=0
        cnam_last='?'
45      read(20+n,'(a80)',end=55)string
        if(string(1:6).ne.'SEQRES') goto 45
        cnam=string(12:12)
        if(cnam.ne.cnam_last) then
          nchain_seqres=nchain_seqres+1
          nres_seqres(nchain_seqres)=0
          cnam_last=cnam
          chain_seqres(nchain_seqres)=cnam
        endif
        goto 45
55      close(20+n)

        nchains_pdb(n)=nchain_seqres
        mychain_pdb(n)=chain_seqres(1)

        if(nchain_seqres.gt.nchain_max) then
          nchain_max=max(nchain_max,nchain_seqres)
          ichain_max=n ! store the ID of the pdb with most chains
        endif

        if(nchain_seqres.gt.1) then
          noligomers=noligomers+1
          write(*,*)'initial_pdb # ',n,' is oligomeric ',nchain_seqres
        endif

      enddo ! do n=1,num_pdbs

C stop if # oligomeric pdb files > 1

      if(noligomers.gt.1) then
        write(*,*)
        write(*,*)'code can only handle *one* oligomeric initial_pdb'
        write(*,*)'it is not yet ready to handle this system'
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

C now, again read the initial_pdb files but this time we get the SEQRES 
C entries - note that we assume that the SEQRES entries of all chains
C are identical (as they should be) so that we only read the first one

      nchain_seqres=1
      nres_seqres(nchain_seqres)=0
      nres_first=0

      do n=1,num_pdbs

C if ialready=1 then we only need to read the SEQRES entries for the
C first chain of the first pdb - there are no other residues to add

        if(ialready.eq.1.and.n.gt.1) exit

        nres_tmp=0
        nres_first(n)=nres_seqres(nchain_seqres)
        open(unit=20+n,file=initial_pdb(n),status='unknown')

50      read(20+n,'(a80)',end=60)string

C skip non-SEQRES lines

        if(string(1:6).ne.'SEQRES') goto 50

C skip SEQRES entries for chains other than the first one

        if(string(12:12).ne.mychain_pdb(n)) goto 50

C otherwise, read and store the residue names...

        do mno=1,13
          jbeg=(mno-1)*4+20
          jend=(mno-1)*4+22
          rnam=string(jbeg:jend)
          if(rnam.eq.'   ') cycle
          nres_seqres(nchain_seqres)=nres_seqres(nchain_seqres)+1
          rnam_seqres(nchain_seqres,nres_seqres(nchain_seqres))=rnam
          nres_tmp=nres_tmp+1
        enddo
        goto 50
60      close(20+n)

        write(*,*)'#res in chain ',n,' = ',nres_tmp

C once we're done with this pdb we'll add on a linker containing
C nlinker_res(:) residues - there's no need to do this for the last 
C pdb file obviously, so we quickly exit if that's the case
     
        if(n.eq.num_pdbs) exit

        do ijk=1,nlinker_res(n)
          nres_seqres(nchain_seqres)=nres_seqres(nchain_seqres)+1
          rnam='   '
          if(nlinker_seq(n)(ijk:ijk).eq.'X') rnam='XXX'
          if(nlinker_seq(n)(ijk:ijk).eq.'A') rnam='ALA'
          if(nlinker_seq(n)(ijk:ijk).eq.'C') rnam='CYS'
          if(nlinker_seq(n)(ijk:ijk).eq.'D') rnam='ASP'
          if(nlinker_seq(n)(ijk:ijk).eq.'E') rnam='GLU'
          if(nlinker_seq(n)(ijk:ijk).eq.'F') rnam='PHE'
          if(nlinker_seq(n)(ijk:ijk).eq.'G') rnam='GLY'
          if(nlinker_seq(n)(ijk:ijk).eq.'H') rnam='HIS'
          if(nlinker_seq(n)(ijk:ijk).eq.'I') rnam='ILE'
          if(nlinker_seq(n)(ijk:ijk).eq.'K') rnam='LYS'
          if(nlinker_seq(n)(ijk:ijk).eq.'L') rnam='LEU'
          if(nlinker_seq(n)(ijk:ijk).eq.'M') rnam='MET'
          if(nlinker_seq(n)(ijk:ijk).eq.'N') rnam='ASN'
          if(nlinker_seq(n)(ijk:ijk).eq.'P') rnam='PRO'
          if(nlinker_seq(n)(ijk:ijk).eq.'Q') rnam='GLN'
          if(nlinker_seq(n)(ijk:ijk).eq.'R') rnam='ARG'
          if(nlinker_seq(n)(ijk:ijk).eq.'S') rnam='SER'
          if(nlinker_seq(n)(ijk:ijk).eq.'T') rnam='THR'
          if(nlinker_seq(n)(ijk:ijk).eq.'V') rnam='VAL'
          if(nlinker_seq(n)(ijk:ijk).eq.'W') rnam='TRP'
          if(nlinker_seq(n)(ijk:ijk).eq.'Y') rnam='TYR'
          rnam_seqres(nchain_seqres,nres_seqres(nchain_seqres))=rnam
        enddo

      enddo ! do n=1,num_pdbs

      write(*,*)
      write(*,*)'after reading SEQRES entries'
      write(*,*)'we find that:'
      write(*,*)
      do n=1,nchain_seqres
        write(*,*)'SEQRES chain ',n,' contains ',
     &            nres_seqres(n),' residue entries'
      enddo

C now write out the complete SEQRES nchain_max times...

      do nnn=1,nchain_max
        do l=1,70
          outstring(l:l)=' '
        enddo
        outstring(1:19)='SEQRES   1 X  123  '
        l=mod(nnn,78)
        if(l.eq.0) l=78
        outstring(12:12)=blah(l:l)      

C make sure SEQRES line has correct #res 

        write(outstring(13:17),'(i5)')nres_seqres(nchain_seqres)

        nline=0
        do mp=1,nres_seqres(nchain_seqres),13
          nline=nline+1
          write(outstring(7:10),'(i4)')nline
          ibeg=mp
          iend=min(nres_seqres(nchain_seqres),mp+12)
          j=0
          do i=ibeg,iend
            j=j+1
            rnam=rnam_seqres(nchain_seqres,i)
            i1=20+(j-1)*4
            i2=22+(j-1)*4
            outstring(i1:i2)=rnam
          enddo

C remember to blank out the rest of the last line

          do j=i2+1,70
            outstring(j:j)=' '
          enddo

          write(12,'(a70)')outstring(1:70)
        enddo
      enddo ! do nnn=1,nchain_max

C finally, read the pdb files again - this time writing out the ATOM and
C HETATM entries but renumbering residues in the ATOM entries - UPDATE
C note that we only do this if ialready.ne.1...

C note that we have to do this nchain_max times, each time incrementing
C the chain ID and each time translating each domain to a different
C location - first, let's figure out how many separate domains we have

      ndomains=num_pdbs*nchain_max
      nonaside=int(real(ndomains)**(1.0/3.0))+1
      disp=5000.0/(real(nonaside))
      icount=0
      jcount=0
      kcount=0

C now we need to know where the center of the (only) oligomer is - we'll
C be displacing all other chains so that they're well away from it...

C --------------------------------------------------------------------
C note that we may want to specifically define the non-moving chain in
C those cases where nchain_max = 1 - at the moment it defaults to the
C first pdb file being the non-mover - which probably should be changed
C using something like the lines written out below...
C --------------------------------------------------------------------

      if(nonmovingpdb.lt.1) then
        if(nchain_max.eq.1) then
          ichain_max=num_pdbs/2+1
        endif
      else
        ichain_max=nonmovingpdb
      endif

      write(*,*)'the non-moving pdb will be # ',ichain_max

C now find the center of geometry for all pdbs...

      do n=1,num_pdbs
        natms=0
        xtot=0.0
        ytot=0.0
        ztot=0.0
        open(unit=20+n,file=initial_pdb(n),status='unknown')
71      read(20+n,'(a66)',end=81)char66
        if(char66(1:4).ne.'ATOM'.and.char66(1:4).ne.'HETA') goto 71
        read(char66,'(30x,3f8.3)')x,y,z
        natms=natms+1
        xtot=xtot+x
        ytot=ytot+y
        ztot=ztot+z
        goto 71
81      close(20+n)
        xcen(n)=xtot/real(natms)
        ycen(n)=ytot/real(natms)
        zcen(n)=ztot/real(natms)
      enddo

C for each pdb file figure out the ID of the first chain - e.g. if the
C chainID is C then the ID will be 3... - this will help (but not
C completely solve) issues with correctly identifying chains

      IDchain_pdb=0
      do n=1,num_pdbs
        do l=1,78
          if(mychain_pdb(n).eq.blah(l:l)) then
            IDchain_pdb(n)=l
            exit
          endif
        enddo
      enddo

C now we want to cycle over nchain_max, repeatedly writing out displaced 
C copies of each pdb file (unless it's the non-mover, in which case we
C write only the chain of interest and don't displace it...)

      ndomains=0
      natoms=0
      do nnn=1,nchain_max

C figure out what the correct chainID we should be looking for in the
C oligomeric pdb file - this should work even if the oligomeric pdb
C starts, say, at chain D *as long as* all other chains are sequential:
C i.e. the order must go D,E,F... not D,F,H etc
C note that since this is only relevant to the oligomeric pdb there is
C no need for a loop over all pdbs here...

        IDchain_req=IDchain_pdb(ichain_max)+nnn-1
        chain_2find=blah(IDchain_req:IDchain_req)

        do n=1,num_pdbs

C if the current pdb is not the one with the max number of chains in 
C it, we increment ndomains and icount,jcount,kcount etc

          if(n.ne.ichain_max) then
            ndomains=ndomains+1
            icount=icount+1
            if(icount.gt.nonaside) then
              jcount=jcount+1
              icount=1
            endif
            if(jcount.gt.nonaside) then
              kcount=kcount+1
              jcount=1
            endif
            if(kcount.gt.nonaside) then
              write(*,*)'error - kcount > nonaside '
              write(*,*)'icount = ',icount
              write(*,*)'jcount = ',jcount
              write(*,*)'kcount = ',kcount
              write(*,*)'nnn    = ',nnn
              write(*,*)'n      = ',n
            endif
          endif

          write(*,*)'check me ',nnn,n,icount,jcount,kcount

          open(unit=20+n,file=initial_pdb(n),status='unknown')

C first, find the first and last residue numbers of the chain we want

          natom=0
65        read(20+n,'(a66)',end=66)char66
          if(char66(1:4).ne.'ATOM'.and.char66(1:4).ne.'HETA') goto 65

C if it's the pdb with the nchain_max number of chains then we need to
C make sure we're reading the right chain at this point
C (for all other pdbs we simply take whatever ATOM/HETA entries we find
C and write them out with the appropriate (new) chainID and with the
C appropriate displaced coordinates)

          if(n.eq.ichain_max) then
            if(char66(22:22).ne.chain_2find) goto 65
          endif
          natom=natom+1
          read(char66,'(22x,i4)')ires
          if(natom.eq.1) ires_frst=ires
          ires_last=ires

          goto 65
66        rewind(20+n)

C now go ahead and read all ATOM/HETA entries for writing out

70        read(20+n,'(a66)',end=80)char66
          if(char66(1:4).ne.'ATOM'.and.char66(1:4).ne.'HETA') goto 70

C if it's the pdb with the nchain_max number of chains then we need to
C make sure we're reading the right chain at this point
C (for all other pdbs we simply take whatever ATOM/HETA entries we find
C and write them out with the appropriate (new) chainID and with the
C appropriate displaced coordinates)

          if(n.eq.ichain_max) then
            if(char66(22:22).ne.chain_2find) goto 70
          endif

C we also should ignore any OXT atoms in if n < num_pbs...
C renumber the atom accordingly
C UPDATE - only do this if it's an ATOM

          if(n.lt.num_pdbs.and.char66(14:16).eq.'OXT'.and.
     &       char66(1:4).eq.'ATOM') then
            write(*,*)'skipping OXT atom in pdb# ',n
            goto 70
          endif

C find the residue 

          read(char66,'(22x,i4)')ires

C if it's the first res or last res of a linker-connected part skip it
C UPDATE - we only do this if the 4th command-line argument = 1

          if(skipfirstlast.eq.1) then
            if(ires.eq.ires_frst.and.n.gt.1) then
              write(*,*)'skipping frst res ',ires,' in pdb# ',n
              goto 70
            endif
            if(ires.eq.ires_last.and.n.lt.num_pdbs) then
              write(*,*)'skipping last res ',ires,' in pdb# ',n
              goto 70
            endif
          endif
 
C renumber the residue accordingly - UPDATE no need to change this line
C for cases where ialready=1 since nres_first will be zero...

          write(char4,'(i4)')ires+nres_first(n)
          char66(23:26)=char4

C otherwise, go on...

          natoms=natoms+1
          write(char7,'(i7)')natoms
          char66(5:11)=char7

C UPDATE - attempted fix to make sure that HETATM entries don't get
C changed to HETA entries accidentally...

          if(char66(1:4).eq.'HETA') char66(5:6)='TM'

C make sure the chainID is set to that of the current chain

          char66(22:22)=blah(nnn:nnn)

C shift each chain by appropriate displacements in x,y,z
 
          read(char66,'(30x,3f8.3)')x,y,z

C if this is not the non-moving molecule then displace coordinates

          if(n.ne.ichain_max) then
            x=x-xcen(n)+xcen(ichain_max)+real(icount)*disp
            y=y-ycen(n)+ycen(ichain_max)+real(jcount)*disp
            z=z-zcen(n)+zcen(ichain_max)+real(kcount)*disp
          endif

          write(char66(31:54),'(3f8.3)')x,y,z
ccc       char66(31:38)=char8

C use the occupancy field to record the "domain" number (i.e. template)

          write(char6,'(f6.2)')real(n)
          char66(55:60)=char6

C set the beta factor field to zero...

          char66(61:66)=' -1.00'
 
C now write out the line

          write(12,'(a66)')char66
          goto 70
80        close(20+n)

        enddo ! do n=1,num_pdbs

C do another chain

      enddo ! do nnn=1,nchain_max

C now, having dealt with all initial_pdb files we'll close final_pdb

      close(12)

      stop
      end


C this is the main code that builds RNC models - it was originally
C conceived as a homology modeling code, so it contains within it the
C possibility of calling sequence alignment codes such as "needle" (part
C of the EMBOSS package) and "promals" (from Nick Grishin's group) -
C these calls are not used by AutoRNC: instead the sequence "alignment"
C is generated in advance using separate code that is called earlier in
C the process - the code also has options to call "dssp", "scwrl", and
C "loopy" but all of those are unnecessary for AutoRNC and so can safely
C be skipped

C please be advised that the code is an *absolute mess* - it certainly
C works, but if you're a hardcore programmer you'll probably be in a
C homicidal rage by the time you've finished wading through it all - the
C messiness is in part due to the above-mentioned change in the purpose
C of the code, and in part due to the author (AHE) being unable to find
C the time to focus on the code 100% - it's the usual story, therefore,
C i.e. "buyer beware" and YMMV etc

C if you run into any major problems running the code please feel free
C to contact me at:

C   adrian-elcock@uiowa.edu


C ---------------------------------------------------------------------

C following are notes that AHE made for himself that log some of the
C major changes that occurred during the development of the code...

C 2023_v1.1

C (1) added support for nfails_chn_max=-1 - use this if we *never* 
C want to extend a chain due to clashes - useful for RNCs where we 
C cannot possibly extend past the last residue in the construct...

C (2) rationalized the nonbonded grid - let's just put all atoms on to
C the grid and make it fit snugly around the template - also let's make
C the cell size be equal to the max of the clash criteria...

C (3) fixed a bug in the clash checking of loopy loops - it seemed to be
C only using clash_intra distance, which is deprecated, so replace with
C the max of the two clash distances...

C (4) use mod-100000 function on "nats" and "nats2" when writing out
C HETATM entries to prevent them from getting ***** for atom #s - see
C format statement 715


C new version for 2022; this is based on:
C for new changes look for:

C 2022 v1.0

C one change that needs documenting and rethinking fully is the
C replacement of CYS SG atoms in scwrl_out with those from scwrl_in - we
C might want to make that a more generally exploitable feature

C UPDATE - probably the best way to proceed is to allow scwrl to be
C skipped entirely - we don't need it for the autoRNC work anyway...
C so introduce new logical:

C skip_scwrl     - ensures that we don't use scwrl
C skip_dssp      - ensures that we don't use dssp 
C skip_loopy     - ensures that we don't use loopy

C this should be used with the following logical:

C skip_all_loops - ensures that we don't use loopy



C ahemodel_2018_v4.7.f

C here we will probably leave as much of the code unchanged as possible
C but we will add paths for all codes that are definitely needed and we
C will remove the (unnecessary) path for .ss2 files: these should be in
C the same folder as the .fasta files and the TEMPLATE_ahemodel.pdb

C we need paths directly to the executables for:
C                                                              checked?

C dssp                  - dssp_path      len_dssp    ! OBLIGATORY YES
C    /home/LAB/BIN/dssp-2.0.4-linux-amd64

C emboss                - emboss_path    len_emboss  ! OPTIONAL   YES
C    /home/LAB/BIN/EMBOSS-6.6.0/needle

C promals               - promals_path   len_promals ! OPTIONAL
C    /home/LAB/BIN/promals_package/bin/promals 
C    /home/LAB/BIN/promals_uni2019/bin/promals

C scwrl                 - scwrl_path     len_scwrl   ! OBLIGATORY YES
C    /home/LAB/scwrl4/Scwrl4

C loopy                 - loopy_path     len_loopy   ! OBLIGATORY 
C    /home/LAB/loopy/test/loopy

C various AHE utilities - code_path      len_code    ! OBLIGATORY

C definitely remove:

C use_SS2_constraints
C truncate_template

C we also need to figure out our final approaches to:

C to fix issue with scwrl f***ing up our tripalmitoylated cysteine we
C are going to have SG atoms in scwrl_out.pdb overwritten with those
C from scwrl_inp.pdb - i.e. retain the original coordinates...

C ss2 files - keep them - we need them for helix/sheets in nascents
C rebuild_tails
C truncate_tails
C promals - keep it - probably useful for us to retain


C new version v4.7:

C (1) makes sure that all atom numbers of fake Glys are <100000
C     this was causing the find_stretch code to screw up
C (2) adds new option on first input line - sidechain_mode
C     this integer controls how dipeptide sidechain atoms are treated
C     =0 ! use all sidechain atoms
C     =1 ! use only CB atoms
C     =2 ! use no sidechain atoms
C (3) adds timers so we know how long it takes to make all models
C    (see ktime1 & ktime2 & ru)
C (4) LATE FIX!!!! discovered that find_clashes_in_loopy_pdb... was
C     being called with clash_distance=0.0 - this was because we have
C     replaced clash_distance with separate clash_intra and clash_inter
C     terms - now fixed this bug - see v4.7 BUGFIX
C (5) LATE FIX!!!! not really a bug but we should also fail if nlink=0
C     when we check all CA-CA distances in the loopy pdb file...
C     see v4.7 BUGFIX

C new version v4.6:

C (1) remove "pause" line from superposition
C (2) add option of putting all HETATMs in a separate file
C (3) use_norestrict allows unrestricted dipeptide confs instead
C     of both_coils...
C (4) made sure that when HETATMs move with rebuilt domains that 
C     they get the correct chainIDs:

C i.e. changed:
C    &                      cur_chainID(my_frag), ! v4.6 make sure
C to:
C    &                      fin_chainID(my_frag), ! chainID is final

C Unfortunately this now also makes clear that such HETATMs will be
C written out in the order that they were added - could be chain B
C first, then D then A then C - this will change every time a new model
C is generated - which will probably throw VMD for a loop - so we need a
C better solution eventually...

C AT LINE 6717 WE NEED TO REORDER HETATMS SO THEY'RE BACK2ORIGINAL
C ALSO LOOK FOR "jx," AND FIX LINES WHERE WE ARE ERRONEOUSLY RENUMBERING
C THE RESIDUES FOR HETATM ENTRIES - FINE TO DO THIS FOR ATOMS BUT SILLY
C TO DO THIS FOR HETATMS  
C UPDATE THE jx, BIT SHOULD BE DONE NOW...


C UPDATE - THIS NEARLY WORKS - THERE'S STILL A LOT TO FIGURE OUT E.G.

C BUT A HUGE ISSUE WITH SYSTEM CALLS TO COPY COMMANDS - THEY SEEM TO
C SCREW UP IN THE TAIL BUILDING PART - NOT SURE WHY BUT I'VE LOST
C CONFIDENCE IN THE SYSTEM CALLS SO I'LL BE ELIMINATING THEM WHEREVER
C POSSIBLE - THIS IS ALREADY DONE IN THE TAIL BIT BUT NOW NEED TO DO
C THIS ALSO IN THE LOOP BUILDING PART...

C see this test directory:

C /home/aelcock/2017_ECOLI_CELL/POLYSOME_MODELING_2020/FRAGMENT_TEST_2PUB_2020

C (1) WHY IS THE FINAL RESIDUE BEING WRITTEN OUT IN THE REMARKS
C     SECTION OF THE FINAL AHEMODELPDB
C     update- FIXED!!!

C (2) WE NEED TO DELETE FAKE GLYS IN LOOPY-INP THAT HAVE BEEN REPLACED
C BY THE TAIL=ASSOCIATED HETATMS
C     this may not be an issue as they are miles away from anything

C (3) WE NEED TO FIGURE OUT HOW THIS WOULD WORK WITH MULTIPLE
C DOMAINS HAVING MULTIPLE SETS OF HETATMS

C (4) WE NEED TO MAKE WORK WITH BACKTRACKING...



C new version v4.5:

C (1) fixes issue with scwrl reordering atoms by latin alphabet in
C sidechain - this is so silly it's unbelievable
C (2) will add the HETATM w/ domain feature
C (3) comment out some unnecessary writes - look for cv4.5

C new version v4.4:

C (1) revamps the nbacktrack mess
C     works *much* better
C     "iworkoutwards" - new logical if "yes" then we try to remove
C     structure from the to-be-built tail
C (2) adds possibility of relaxing clash criteria for first few
C     residues of a tail - ONLY IMPLEMENTED FOR NTER TAILS...
C     "nrelax" - # of residues to relax criteria for
C     "frelax" - factor by which to scale clash distances

C new version v4.3:

C (1) adds nbacktrack - # of residues to backtrack when a failure to
C     build a tail occurs...
C     set nbacktrack=0 for backwards compatibility with earlier versions

C UPDATE - THIS IS DEFINITELY NOT FINISHED BUT IT SORT-OF WORKS...

C new version v4.2:

C (still don't have a solution to the problem of directly associating
C HETATM entries with moving protein chunks)

C (1) add different clash cutoffs for ATOM-ATOM & ATOM-HETATM clashes
C     change clash_distance in input to clash_dist1, clash_dist2

C where clash_dist1 is for ATOM-ATOM --> clash_aa2
C       clash_dist2 is for ATOM-HETATM --> clash_ha2

C adds in ihet_grid = 1 for HETATMs on grid =0 for ATOMs on grid

C (2) puts in a fix for when the .ss2 file is longer than the .fasta -
C we don't want this to be a fatal error...

C new version v4.1:

C (1) we will add option for call to promals to use updated, i.e. large
C sequence databases

C this means using a new logical:

C       use_promas3D_2019

C which then calls code in the following directory:

C       /home/LAB/BIN/promals_uni2019/

C instead of:

C       /home/LAB/BIN/promals_package/

C (2) put in a fix so that non-standard entries in uniprot fasta files
C are converted to standard form - this is because of entries like
C P24183 which contains a "U" residue

C (3) add in a fix for lipidation proteins so that we don't
C automatically increase the size of loops when there's a single
C N-terminal or C-terminal residue...

C (4) add in code so that HETATM entries are only added to the
C scwrl_frame.pdb if they are within 6A of an ATOM entry - this is
C because scwrl seems to cock up if the #atoms in the frame file is
C >32500 or so - I found this when trying to run AHEMODEL on 5UYM - this
C seems really crap but I don't see that we can do anything about it...
C one possibility would be to sort all HETATMs according to their
C closest distance and then keeping only the closest 32,000 or so - that
C would probably work but I don't really need it right now...

C (3) we will finish implementation of HETATM entries so that they can
C be placed during building of tails

C (4) we will convert the code into a useable form for outside groups 


C new version v4.0: do the following

C (1) change to consider adding is to add code so that it won't quit if
C a sequence doesn't match to a pdb - we could make the code find a
C match for each sequence even if it's not the absolute best that could
C be found for that particular chain - this could help us deal with ABC
C transporters where we often need to match two sequences to the same
C type of chain...

C UPDATE: DONE

C (2) set up automatic calling of MEMOIR for membrane chains - need to 
C read a file that states localization of each uniprot_ID
C for this, we will introduce align_mode = 4

C UPDATE: NOT ATTEMPTED AND WILL NOT BE ATTEMPTED

C new version v3.9: DELETES CHANGES MADE IN v3.8

C introduces new logical : have_ss2_file

C new version v3.8: DEAD

C new version v3.7:

C puts HETATM entries into a file for scwrl and uses them with
C -f option...

C new version v3.6:

C only moves the promals3D folder if there is a different fasta to run
C (not done yet)
C allows loops to not be built - skip_all_loops (done)
C fixes bug with submissions to promals3D if the input fasta doesn't
C contain a return character 
C no need now for .ss2 file to be in a subdirectory

C new version v3.5:

C implements

C sheet_cut +

C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Nter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Cter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_sheet.pdb

C also allows id_thr to be set for promals3D


C new version v3.4:

C implements 

C helix_cut +

C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Nter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Cter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_helix.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_coils.pdb


C new version v3.3:

C same as parent but removes reading of flatfile - we're only using this
C to get the PSIPRED prediction but we can get this from .ss2 anyway...

C new version v3.2:

C same as parent but fixes issues with array bounds being exceeded

C --------------------------------------------------------
C need to fix array sizes when residue counts exceed 50000
C --------------------------------------------------------

C new version v3.1:

C is currently set to run with the athomas dipeptide.pdb files

C new version v3.0:

C is currently set to run with the PISCES dipeptide.pdb files
C add a path for the ss2 file - this means changing the input file but
C allows us to use different ss2 files - remember that a uniprot of
C PXXXXX means no looking for a ss2 file
C ss2_path

C reads nmodels_desired from input_file

C selects helical residues according to PSIPRED prediction confidence
C seems to work well


C new version v2.9:

C this is a continuation of v2.7...

C (1) add in write-out of PSIPRED confidence to CARTOON.pdb instead of
C domain IDs

C note that our method for building N-terminal and C-terminal tails
C appears similar to that described in:

C Efficient conformational ensemble generation of protein-bound peptides
C by Yan, Zhang, and Huang - JCheminform 2017 v9 p59

C new version v2.7:

C (1) avoid writing growing tails to pdb - DONE
C (2) put tail clashes on grid
C (3) allow loop to make multiple models - DONE
C     this requires that we store loop and tail defos so that 
C     they can be restored to their original values...
C (4) write out alignment stats - DONE

C new version v2.6:

C randomly order the N-terminal and C-terminal tails...

C Here is the issue with loopy (I think) - we're making a smaller
C loopy_inp.pdb to speed the calculations up, but I bet we're not then
C recombining the loopy-generated tail with the original loopy_inp.pdb
C to make a complete pdb file - it's easy to fix but...
C Actually that seems like bullshit - perhaps it's more to do with a
C failure in loopy not being recognized for what it is - probably best
C to make sure that we keep a copy of every successful loop-built
C structure so that we can debug line by line...

C goals for Wednesday:

C (1) streamline generation of temp4_XXXX_XXXX files - do this once only
C for N-ter and C-ter and use chainID to distinguish redundancies
C (2) read combined dipeptide files

C goals for Monday:

C (1) reconfigure chainID and res numbering throughout
C (2) implement bounding-box truncated loopy
C (3) implement nonbonded-grid tail building
C (4) clean-up any issues with C-terminal OXT atoms...
C     simple workaround is to delete the final residue!

C new version v2.1:

C various bugfixes required for the use of multiple chains when ires
C exceeds 10,000 - as it typically will do - it's all doable, just a
C bit of a pain to implement.. - in particular, we use tot_res counts a
C lot but this fails once we get past 10,000 - on the other hand, if we
C do it all by chain-local numbering we need to worry about cases where
C we have more than 78 chains - need a solution to this now!

C we need a complete overhaul :(

C scwrl can be given recurring chainIDs with chain-local res# - it
C ignore residues in a chain if they pass from 9999 to 0

C UPDATE we probably should do all of these using each chainID as much
C as possible before moving on to the next: e.g. A 1-100, A101-250, etc
C this is because we have definitively shown that scwrl will reorder
C atoms by chain - so if the same chainID appears multiple times in
C scwrl_inp.pdb the atoms with the same chainID will be lumped together
C in scwrl_out.pdb. 

C SO THE BOTTOM LINE APPEARS TO BE THAT WE SHOULD USE THE SAME APPROACH
C FOR SCWRL, LOOPY, AND TAILS - NAMELY - USE EACH CHAINID AS MUCH AS
C POSSIBLE BEFORE MOVING ON TO THE NEXT ONE - I.E. USE CHAIN A
C CONTINUALLY BUT WITH INCREASING RESIDUE NUMBERS UNTIL 9999 IS GOING TO
C BE EXCEEDED BY ONE OF THE CHAINS - IF SO, THEN THAT CHAIN SHOULD BE
C GIVEN THE NEXT CHAINID ETC...

C --------------------------------------------------------------
C THERE IS A LOT TO FIGURE OUT HERE BUT I THINK IT WILL ALL WORK
C --------------------------------------------------------------

C loopy must be given recurring chainIDs with chainID-specific res#
C tails must be given recurring chainIDs with chainID-specific res#

C UPDATE - IT SEEMS CERTAIN THAT LOOPY WILL GO MUCH FASTER IF ONLY GIVEN
C RELEVANT ATOMS - IT SEEMS LIKE IT NEEDS THE BACKBONE ATOMS OF RESIDUES
C TO BE DEFINED SO SAFEST TO BUILD A BOUNDING BOX AROUND THE BEG AND END
C OF THE LOOP AND THEN KEEP ENTIRE RESIDUES IF THEY HAVE ANY ATOM WITHIN
C THE BOUNDING BOX - VERY SIMPLE TO IMPLEMENT...
 
C new version v2.0:

C (1) bugfix for i_am_original in Nter & Cter tails
C (2) rebuild_auto
C (3) use iam_in_pdb_chn to store domain #
C     note that the default value is now -1 - this is for missing
C     residues - residues in loops will get domain # 0...
C (4) chainID_desired


C need to do the following:

C (1) make sure that loop merging/expansion etc is applied only once to
C     each type of chain (assuming this is possible structurally) - at
C     the moment different copies of the same chain can expand in
C     different ways, which makes no sense
C     DONE DONE
C
C     NOTE THAT THE ABOVE MAY NEED TO BE MADE OPTIONAL - I THINK IT
C     SHOULD WORK IN ALL CASES AS WE ONLY APPLY IT IN CASES WHERE THE
C     INITIAL LOOP ASSIGNMENTS ARE IDENTICAL, BUT I WORRY THAT THERE
C     MIGHT BE PROBLEMS IN CASES WHERE THE CHAINS DIFFER IN TERMS OF
C     MISSING RESIDUES - THIS IS POSSIBLY AN UNFOUNDED WORRY...
C
C (2) need to apply the SS constraint to all parts of loop merging -
C     it's not currently applied in the very earliest stages
C     DONE DONE
C
C (3) cosmetically, need to add the loop sequence when writing out the
C     FINAL FINAL LOOP statements
C     DONE DONE
C 
C (4) need to allow very long loop regions to be skipped entirely - this
C     would greatly simplify the use of a TEMPLATE that in reality
C     contains more than one (and physically dispersed) templates...
C
C (5) need separate code to merge TEMPLATE_ahemodel.pdb files...

C new version v1.9:

C same as parent (v1.8) but:

C (1) attempts to merge loops immediately (does again later still)
C (2) only writes "distance check" for crazy bonds

C new version v1.8:

C same as parent (v1.7) but:

C (1) allows ss_weight to be set for promals3D
C (2) does clash checking on loops
C (3) allows tails to not be built
C     build_no_tails
C (4) when extending loops gives preference to SS elems

C     integer i_in_SS_elemnt_tot(10000) ! for all res in all chains

C new version v1.7:

C same as parent (v1.6) but:

C (1) allows use of PROMALS3D for the sequence alignments...

C new version v1.6:

C same as parent (v1.5) but:

C (1) allows us to continue if we find a template chain with no match -
C     in previous versions this is a fatal error - in this version we
C     will allow the possibility of continuing...

C     continue_regardless

C (2) to this end, we also allow resname 'XXX' to be used - this gets
C     convert to 'X' in fastas for needle

C new version v1.5: 

C same as parent (v1.4) but:

C (1) finishes implementation of writing out DSSP assignments
C (2) writes out PSIPRED SS predictions also

C new version v1.5:

C same as parent (v1.4) but:

C (1) when use_SS_constraints="yes" the code first does alignments
C     without SS constraints, then uses the alignments to determine if the
C     template needs to be truncated - this solves problems where the
C     alignment goes crazy with longer-than-necessary templates...
C     note that this truncation only takes effect if we set the logical
C     truncate_template to yes
C     for the code alterations look for v1.4 and also:
C
C       len_of_pdb_chn_fin
C       ifrst_time_thru
C
C (2) writes out stats on %identities etc for each chain - not yet fully
C implemented the gap2 stuff for template residues unaligned...

C syntax for possibly adding path to dipeptide files...

C filename_full='RUN_'//poop//'/'//filename_stem

C new version v1.3:

C same as parent (v1.2) but:

C (1) adds "isentback=1" lines to speed up RMS restart of tails
C     see lines that contain "BUGFIX ADDED LATE TO SPEED UP"
C (2) allows override of ibeg_tot_nters and iend_tot_cters
C     we introduce new logical "truncate_tails" and add entry
C     in input file....
C (3) writes out alignments in rosetta format also
C     generates .rewritten.rosetta files...
C (4) makes comfile for running rosetta *with* symmetry, which needs to
C     edited before running (at the symmdef stage...)

C new version v1.2:

C same as parent (v1.1) but:

C (1) randomly samples from the nloopy confs
C (2) declared membrane_dive as logical (was declared as integer!)


C BUGFIX v1.1
C new bug fix to solve problem when merging loop with Nter or Cter tail
C that contains only one residue...

C this version is identical to its parent but it writes out
C the input file into REMARK entries in the AHEMODEL.pdb

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops_correct_Cter_membrane_aware_hand_align_optional_seqres_mandatory_needle_tail_rebuild_optional_input_file_HETATM.f

C this version is a child of its parent:

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops_correct_Cter_membrane_aware_hand_align_optional_seqres_mandatory_needle_tail_rebuild_optional_input_file.f

C with the following changes:

C (1) HETATM entries are only kept if they are within 4A of any atom in
C the template pdb file - this decision is made in an iterative way...

C (2) occupancy entries are assumed to be zero in the template and if
C they are not then the input number is assumed to be a rigid-domain
C number - we will add to this during the building - occupancy entries
C that are written out are the final rigid-domain numbers

C (3) HETATM entries can be associated with residues if we use
C rebuild_tails...


C this version is a child of its parent:

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops_correct_Cter_membrane_aware_hand_align_optional_seqres_mandatory_needle_tail_rebuild_optional_input_file.f

C with the following changes:

C (1) change beta factor write-out so that 
C     sequence identity + structure = 99.99
C     sequence similarity + structure = 75.00
C     sequence aligned + structure = 50.00
C     no structure = 0.00

C (2) implement an input file instead of command line arguments

C (3) implement a membrane building routine


C this version is a child of its parent:

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops_
C correct_Cter_membrane_aware_hand_align_optional_seqres_mandatory_
C needle.f 

C we will institute the following changes:

C (1) the biggest change will be to allow the option of overriding where
C the N- and C-terminal tails start in a structure - we want to allow
C users to provide multi-domain proteins as spatially separated domains,
C have ahemodel fill in their loops (but not the inter-domain loop) and
C then build all but one of the domains as N- or C-terminal extensions
C of the remaining domain - to switch off building of the inter-domain
C loops we will need to eliminate any loop that contains the over-ridden
C definition of the N- or C-terminal loop - other than that, we'll
C basically just be reusing the dipeptide-building code but, instead of
C getting random dipeptides from the library we'll be getting them from
C the structure that we've just built... 

C STATUS: DONE - we're still confused about how to make sure that
C structured parts are properly assigned beta=99.99 in the final
C structure - we're using i_am_original, i_done_res, inot_made_yet and
C I'm quite confused about what means what - this should get resolved...
C see point (3) below...

C more minor differences to implement:

C (2) when a loop-build fails we should increase the size of the loop by
C only one residue (the least conserved) not by both - DONE

C (3) we would like to retain 0.00 entries (i.e. modeled residues) for 
C residues that were previously built in earlier runs of ahemodel - 
C currently those entries get set back to 99.99 (i.e. structured)  
C because the code doesn't try to track them - this will be important 
C for properly illustrating what's built and not-built for more 
C complicated protein architectures - NOT DONE

C (4) we're at the point where we should get run-time parameters from a
C file instead of from the command line - NOT DONE

C (5) we should add the possibility of associating HETATM entries with
C the nearest bit of structure - this would enable them to be carried
C with them during model construction, so for example we could have the
C bound-DNA used in the construction of a TF...

C this version is a child of its parent:

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops_
C correct_Cter_membrane_aware_hand_align_optional_seqres_mandatory_
C NWalign.f

C that version was an aborted one that used NWalign instead of needle to
C do all the sequence alignments. The advantage of that code was that it
C allowed secondary structure elements found in the template to be used
C as constraints in the Needleman-Wunsch global alignment. The
C disadvantage was that NWalign penalizes end gaps in the same way that
C it penalizes internal (i.e. regular) gaps. This is undesirable as it
C tends to make the template and target line up at the very beginning
C and the very end of the sequence in order to not have end gaps.
C *So*, here we will reinstate the use of needle since: (a) needle does
C not penalize end gaps (in its default operation), and (b) we have
C managed to make a version of needle (embaln.c) that can read in
C secondary structure elements (from SS_elements.txt) and use these as
C constraints in the global alignment...

C this version is a child of its parent:

C ahemodel_loop6_timeout_deletion_9999res_forbid_unclosed_loops.f

C here we will attempt to add three new features:

C (1) correct handling of OXT atoms that are already correct - original
C     code version deletes these for some reason - DONE!!!
C (2) reading of membrane coordinates so that tails and loops are aware
C     of the membrane and will try to avoid it - this may require a
C     command-line argument so that tails/loops can be allowed to enter
C (3) option to use hand-made needle alignments and to stop immediately
C     after making the alignments - needs command line argument - DONE!!!
C (4) need to add extra loops where the alignment *thinks* it's great
C     but it actually has a massive distance between two adjacent
C     residues in the target sequence - this can happen when the
C     sequence similarity is poor and there is a big gap in the template
C     with the residues either side of the gap both matching adjacent
C     residues in the target sequence
C (5) allow option to get fasta directly from SEQRES entries in pdb

C unresolved issues:

C (1) possibly add an initial pdb check to eliminate residues whose
C     backbones are not fully resolved - this may not be a real problem
C     - it certainly *appeared* to be an issue for scwrl, but this later
C     turned out to be due to alternate A/B conformers in the original
C     pdb file screwing things up - the code has been updated so that it
C     ignores any B or other conformers and so far we have not had
C     trouble with incomplete backbone atoms in any pdbs - there's
C     always a first time however...
C (2) need to make sure the code will work with num_res>9999
C (3) if making loops bigger (due to clashes) then we may need to ensure
C     that we properly merge with other loops or tails - this may also,
C     however, not be necessary
C (4) need to add in OXT atoms to the final model
C (5) quit after reading first model in case of NMR structures...

      program place_molecules
c
      implicit real(a-h,o-z)
      parameter   (pi=3.141592654)

C v2.1 set this as a parameter for testing purposes only - final value
C will be 78 probably - possibly less if some are illegal characters

      integer, parameter :: chainIDmax=78

C v3.2 set max# of residues as a parameter - set to 1 million...

      integer, parameter :: ntotresmax=1000000

      INTEGER, PARAMETER :: K4B=selected_int_kind(9)

      interface
      subroutine pdbsup(nat_loc,xm_loc,ym_loc,zm_loc,rf_loc,rms_loc)

      integer, intent(in)     :: nat_loc
      real, intent(inout)     :: xm_loc(1:nat_loc)
      real, intent(inout)     :: ym_loc(1:nat_loc)
      real, intent(inout)     :: zm_loc(1:nat_loc)
      real, intent(in)        :: rf_loc(1:4,1:3)
      real,    intent(out)    :: rms_loc

      end subroutine
      end interface


      integer      :: omp_get_thread_num

C v3.4 store first and last conformation for each type of SS dipep

      integer i_need_dipep(20,20)
      integer dipep_both_coils_beg(20,20)
      integer dipep_both_coils_end(20,20)
      integer dipep_helix_Nter_beg(20,20)
      integer dipep_helix_Nter_end(20,20)
      integer dipep_both_helix_beg(20,20)
      integer dipep_both_helix_end(20,20)
      integer dipep_helix_Cter_beg(20,20)
      integer dipep_helix_Cter_end(20,20)
      integer dipep_sheet_Nter_beg(20,20)
      integer dipep_sheet_Nter_end(20,20)
      integer dipep_both_sheet_beg(20,20)
      integer dipep_both_sheet_end(20,20)
      integer dipep_sheet_Cter_beg(20,20)
      integer dipep_sheet_Cter_end(20,20)
    
C v4.5 
     
      character*80 char80
      character*30 char30
      character*66 char66
      character*24 char24

      real, allocatable :: xcys(:)
      real, allocatable :: ycys(:)
      real, allocatable :: zcys(:)

      character*7 greek ! stores first 7 letters of greek alphabet
 
      character*1 ajunk9
      character*1 atail
      character*1 chainID
      character*1 chainID_final
      character*1 chainID_desired
      character*1 chainID_current
      character*1 chainID_last
      character*4 atmnam   
      character*4 char4
      character*20 char20
      character*46 tailstring
      character*26 tailname
      character*27 loopname
      character*22 loopnam2
      character*6 char6
      character*9 char9
      character*70 outstring
      character*3 rnam 
      character*3 resnam 
      character*89 string_sup
      integer     resnum,num_box

C v2.7
C string66_tail is the line that will be written to the tail pdb file
C eventually - up to 1000000 atoms allowed
C ires_tot_tail is this tail atom's total residue number...

      character*30              :: string30
      character*66              :: string66_tail(1000000)
      character*66              :: string66_tail_tmp(1000000)
      integer                   :: ires_tot_tail(1000000)
      integer                   :: ires_tot_tail_tmp(1000000)
      integer                   :: dtail(1000000) ! v4.5
      real                      :: xtail(1000000)
      real                      :: ytail(1000000)
      real                      :: ztail(1000000)
      integer                   :: dtail_tmp(1000000) ! v4.5
      real                      :: xtail_tmp(1000000)
      real                      :: ytail_tmp(1000000)
      real                      :: ztail_tmp(1000000)

C v3.0 

      character*3               :: rnam1,rnam2
      character*1               :: char1
      character*9999            :: dssp_tmp
      character*1000 grepstring

C v2.3

      character*3   char3,char3b,char3c,char3d
      character*1000 dipepname
      character*66, allocatable :: string66_dipep(:,:,:,:)
      real,         allocatable :: x_dipep(:,:,:,:)
      real,         allocatable :: y_dipep(:,:,:,:)
      real,         allocatable :: z_dipep(:,:,:,:)

      character*66         :: string66_first
      character*66, allocatable :: string66_temp4(:,:)
      integer,      allocatable :: n_temp4(:) ! # atoms in res
      integer,      allocatable :: d_temp4(:,:) ! v4.5 domain#
      real,         allocatable :: x_temp4(:,:)
      real,         allocatable :: y_temp4(:,:)
      real,         allocatable :: z_temp4(:,:)
      integer                   :: temp4_N(9999)
      integer                   :: temp4_CA(9999)
      integer                   :: temp4_C(9999)
      integer                   :: temp4_O(9999)

      real              :: rf(4,3)
      real              :: rf_prev(4,3)
      integer dipep_natoms(20,20)
      integer dipep_natom1(20,20) ! atoms in first res
      integer dipep_natom2(20,20)
      integer dipep_N1(20,20)
      integer dipep_N2(20,20)
      integer dipep_CA1(20,20)
      integer dipep_CA2(20,20)
      integer dipep_C1(20,20)
      integer dipep_C2(20,20)
      integer dipep_O1(20,20)
      integer dipep_O2(20,20)

C v4.6 

      integer ichain(1000000)

C 4.7

      integer*8   ktime1 
      integer*8   ktime2
      integer*8   ru

C v2023_v1.1

      integer*8   nfails_chn

C v2.1

      integer     num_chainID
      integer     nres_this_chainID
      character*1 cur_chainID(10000)
      character*1 fin_chainID(10000)
      integer     cur_offset(10000)
      integer     fin_offset(10000)
      real        xq(2),yq(2),zq(2)
      character*1 cwithin(100000)
      integer     rwithin(100000)
      integer     lwithin(100000)
      integer     chainlocal_2_totres(chainIDmax,10000)

      integer, allocatable     :: temp_grid(:,:,:)   ! stores #in box
      character*1, allocatable :: cchn_grid(:,:,:,:) ! cur_chainID
      integer, allocatable     :: rchn_grid(:,:,:,:) ! res#
      integer, allocatable     :: ihet_grid(:,:,:,:) ! =1 if HETATM
      real,    allocatable     :: xatm_grid(:,:,:,:) ! xcoord
      real,    allocatable     :: yatm_grid(:,:,:,:) ! ycoord
      real,    allocatable     :: zatm_grid(:,:,:,:) ! zcoord
  
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
c     real, allocatable :: rf(:,:)
c     real, allocatable :: rm(:,:)

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
      character*6              :: chain_nam(1000)
      integer                  :: nres_totl(1000)
      integer                  :: nres_iden(1000)
      integer                  :: nres_posi(1000)
      integer                  :: nres_alin(1000)
      integer                  :: nres_gaps(1000)
      integer                  :: nres_befr(1000)
      integer                  :: nres_aftr(1000)
      integer                  :: nres_tmplte_gaps(1000)
      integer                  :: nres_tmplte_befr(1000)
      integer                  :: nres_tmplte_aftr(1000)
      character*50000,allocatable :: seq_of_fasta(:)
      character*50000             :: seqss2 ! v3.3
      character*60                :: seq_chunk
      character*50000,allocatable :: seq_of_pdb_chn(:,:)
      integer,        allocatable :: iam_in_pdb_chn(:,:,:)
      integer,        allocatable :: iam_in_SS_elem(:,:,:)
      character*1,    allocatable :: my_DSSP_design(:,:,:)
      character*1,    allocatable :: my_PSIPRED_prd(:,:,:)
      character*1,    allocatable :: my_PSIPRED_cnf(:,:,:) ! confidence
      real,           allocatable :: my_PSIPRED_vl1(:,:,:) ! coil conf
      real,           allocatable :: my_PSIPRED_vl2(:,:,:) ! helx conf
      real,           allocatable :: my_PSIPRED_vl3(:,:,:) ! shet conf
      integer        ,allocatable :: len_of_pdb_chn(:,:)
      integer        ,allocatable :: len_of_pdb_chn_fin(:,:)
      integer        ,allocatable :: num_res_in_seq(:)
      character*1                 :: chn_last
      character*22                :: loopy_file
      character*100               :: string100

      integer num_domains(100000) ! number of rigid domains in frag
      integer num_res_in_frag(100000)
      integer my_frags_seq_num(100000)
      integer idone(50000)
      integer ibeg_word(50000)
      integer iend_word(50000)
      integer ibeg_loop(50000)
      integer iend_loop(50000)
      integer ibeg_ori_loop(50000)
      integer iend_ori_loop(50000)
      integer ibeg_loop_final(50000)
      integer iend_loop_final(50000)
      integer i_in_pdb_fasta(50000) ! only for res in current chain
      integer i_in_SS_elemnt(50000) ! only for res in current chain
      integer i_in_pdb_coord(50000)
      integer i_am_in_a_del(50000)
      integer my_nnn(50000)
      integer my_mmm(50000)
      integer my_nres2(50000)
      integer my_nres2_tru(50000)

C v4.7
 
      integer sidechain_mode

C v1.2

      integer imad(10)
      integer jmad(10)
      integer flen


      integer num_tot_res
      integer num_tot_nters
      integer num_tot_cters
      integer num_tot_loops
      integer my_frst_res(100000) ! dimensioned 1:# total chains
      integer my_last_res(100000)
      integer my_totl_res(100000)
      integer irebuild_nters(100000)
      integer irebuild_cters(100000)

      character*1 ichn_tot_nters(100000)
      integer ifrg_tot_nters(100000)
      integer ibeg_tot_nters(100000)
      integer iend_tot_nters(100000)
      integer ibeg_loc_nters(100000)
      integer iend_loc_nters(100000)
      character*1 ichn_tot_nters_original(100000)
      integer ifrg_tot_nters_original(100000)
      integer ibeg_tot_nters_original(100000)
      integer iend_tot_nters_original(100000)
      integer ibeg_loc_nters_original(100000)
      integer iend_loc_nters_original(100000)

      character*1 ichn_tot_cters(100000)
      integer ifrg_tot_cters(100000)
      integer ibeg_tot_cters(100000)
      integer iend_tot_cters(100000)
      integer ibeg_loc_cters(100000)
      integer iend_loc_cters(100000)
      character*1 ichn_tot_cters_original(100000)
      integer ifrg_tot_cters_original(100000)
      integer ibeg_tot_cters_original(100000)
      integer iend_tot_cters_original(100000)
      integer ibeg_loc_cters_original(100000)
      integer iend_loc_cters_original(100000)

      character*1 ichn_tot_loops(100000)
      integer ifrg_tot_loops(100000)
      integer ibeg_tot_loops(100000)
      integer iend_tot_loops(100000)
      integer ibeg_loc_loops(100000)
      integer iend_loc_loops(100000)
      integer ibeg_ori_loops(100000)
      integer iend_ori_loops(100000)

      integer num_act_cters          ! # of actual already-placed Cters
      real rand_res(1000000)

C v2.6:

      integer typ_tot_tails(100000)
      integer imadtail(100000)
      integer jmadtail(100000)

      real a(3),b(3),c(3),vca(3),vcb(3),xp(3),yp(3),zp(3)

C v3.9:

      logical have_ss2_file
      logical ss2_file_exists

C v4.4:
      logical iworkoutwards

      logical skip_all_loops
      logical use_SS_constraints
      logical rebuild_tails
      logical rebuild_auto
      logical build_no_tails
      logical truncate_tails
      logical truncate_template ! only used if use_SS_constraints true
      logical rosetta_symmetry
      logical continue_regardless
      logical use_promals3D ! this overrides use_SS... & truncate_temp.
      logical use_promals3D_2019 ! uses uniref90_2019 version 
      logical use_tmalign   ! set to "no" to for no 3D info in promals

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
      real, allocatable :: xs(:) ! atoms in scwrl_out.pdb
      real, allocatable :: ys(:)
      real, allocatable :: zs(:)
      character*66, allocatable :: sh(:) ! v4.5
      real, allocatable :: xh(:) ! HETATM entries already there
      real, allocatable :: yh(:)
      real, allocatable :: zh(:)
      character*66, allocatable :: si(:) ! v4.5
      real, allocatable :: xi(:) ! HETATM entries already there
      real, allocatable :: yi(:)
      real, allocatable :: zi(:)

C v4.5 9 April 2021 - we keep a copy of the ih array (ih_original) which
C we use to restore ih after we've completed each MODEL - if we don't do
C this then we forget to look for attached HETATM entries when we build
C in structured domains as part of tails...

      integer, allocatable :: ih(:) ! keep this HETATM entry
      integer, allocatable :: ih_original(:) ! keep this HETATM entry
      integer, allocatable :: mh(:) ! original HETATM#
      integer, allocatable :: mi(:) ! original HETATM#
      integer dm(10000) ! v4.5 domain# of atom to be added
      real xm(10000) ! amino acid to be added
      real ym(10000)
      real zm(10000)
      integer iresm(10000)
      real xn(10000) ! amino acid to be added
      real yn(10000)
      real zn(10000)
      integer iresn(10000)
      real xc(10000) ! amino acid adjacent to one to be added
      real yc(10000)
      real zc(10000)
      real xl(10000) ! amino acids close to the one to be added
      real yl(10000)
      real zl(10000)
      real xll(10000) ! used for checking whether loopy loops are closed
      real yll(10000)
      real zll(10000)
      real xlast(1000000)
      real ylast(1000000)
      real zlast(1000000)

C v3.2 make many arrays dependent on the parameter ntotresmax

      character(len=ntotresmax) seq_tot_res
      character(len=ntotresmax) chn_tot_res
      character(len=ntotresmax) seq_scwrl
      integer res_tot_res(ntotresmax)
      integer ipos_act_cters(ntotresmax) ! scwrl res# CTE
      integer my_scwrl_res(ntotresmax)
      integer my_tot_res(ntotresmax)
      integer i_done_res(ntotresmax)
      real    my_bfactr(ntotresmax)   ! per-res based on i_am_original
      integer my_domain(ntotresmax) ! per-res domain #...
      integer i_old_unstruc(ntotresmax) ! =1 if unstruc in original pdb
      integer i_am_original(ntotresmax) ! =1 if present in original pdb
      integer i_am_original_original(ntotresmax) ! =1 if present in original pdb
      integer i_am_by_loopy(ntotresmax) ! =1 if made by loopy
      real    occ_psipred(ntotresmax)
      integer i_in_pdb_fasta_tot(ntotresmax) ! for all res in all chains
      integer i_in_dom_fasta_tot(ntotresmax) ! for all res in all chains
      integer i_in_SS_elemnt_tot(ntotresmax) ! for all res in all chains
      integer inot_made_yet(ntotresmax)
      integer inot_made_yet_original(ntotresmax)



      character*12 my_fmt
      character*80 my_fmt2
      character*20 my_fmt3
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
      character*10000 string
      character*15    temp4_nter
      character*15    temp4_cter
C     character*100   fasta_get ! v4.5
      character*200   fasta_get
      character*100   wcount
      character*100   touch 
      character*66    betastring
      character*47    sys_stringw
      character*35    sys_stringy
      character*77    sys_stringx
      character*1000  loopstring  ! v1.8 loopy
      character*1000  pro_string0 ! v1.7 promals3D
      character*1000  pro_string1 ! v1.7 promals3D
      character*1000  pro_string2 ! v1.7 promals3D
      character*1000  pro_string3 ! v1.7 promals3D
      character*1000  pro_string4 ! v1.7 promals3D
      character*1000  sys_string1
      character*1000  sys_string2
      character*1000  sys_string3
      character*1000  sys_string4
      character*1000  sys_string5
      character*1000  sys_string6
      character*1000  sys_string7
      character*1000  sys_string8
      character*1000  sys_string9
      character*1000  sys_string9b ! used for VMD cartoon
      character*1000  sys_string9c ! v4.6 used for fixed HETATMs
      character*12    fasta_tmp
      character*80    shit1
      character*80    shit2
      character*80    shit3
      character*80    shit4
      character*80    strong
      character*80    strang(100000)
      character*80    strrng(100000)
      character*80    last_strong
      character*100 strung1
      character*100 strung2
      character*100 strung3
      character*100 strung4
      character*100 strung5
      character*90 streng
      character*90 streng0
      character*90 streng1
      character*90 streng2
      character*4  num_of_res
      character*4  anam
      character*1  cc
      character*1  b1,b2,b3
      character*75 a1,a2,a3,a4,a5,a6,a7,a8,a9
      character*75             :: falign_tmp
      character*75,allocatable :: falign(:,:,:)
      character*75,allocatable :: falign2(:,:,:) ! for needle rewrite
      integer,     allocatable :: ialign(:,:,:)
      real,        allocatable :: salign(:,:,:)
      real,        allocatable :: salign1(:,:,:) ! used when 
      real,        allocatable :: salign2(:,:,:) ! truncate_template
      real,        allocatable :: salign3(:,:,:) ! = 'yes'
      character*78 blah
      character*20 dipep
      character*2  dipep_cur
      character*3  nam_het_typs(100000)
      integer      num_hets_found
      character*3  nam_hets_found(100000)
      integer num_peps(1:20,1:20)
      integer align_mode
      logical membrane_dive
      logical separate_hetatms
      logical use_norestrict
      logical skip_scwrl
      logical skip_dssp 
      logical skip_loopy

      character*1 last_ss_type
      character*1  curr_ss_type
      character*1  my_dssp_type(8)
      character*1  my_res_dssp_type(1000000) ! dssp type of each res
      integer      nss(1000000)
      integer      do_me_still(1000000)
      integer      ibeg_ss(1000000)
      integer      iend_ss(0:1000000)
      integer      ilen_ss(1000000) ! length of SS element
      integer      ibeg_ss_tmp(1000000)
      integer      iend_ss_tmp(1000000)
 
      integer      tmplte_restyp(1000000)
      integer      target_restyp(1000000)
   
      real         id_thr ! used for promals3D
 
      character*1          :: sequenceA(100000)
      character*1          :: sequenceB(100000)
      character*1          :: sequenceM(100000)
      character*1          :: sequenceS(100000) ! SS element
      character*1          :: sequenceD(100000) ! DSSP designation
      character*1          :: sequenceP(100000) ! PSIPRED prediction
      character*1          :: sequenceC(100000) ! PSIPRED confidence
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

C v3.1

      integer, allocatable :: iwant_this_gly(:)

      character*1000       :: strong1000 
      character*1000       :: string1000 
      character*1000       :: dipeptide_path
      character*1000       :: ss2_path
      character*1000       :: dssp_path
      character*1000       :: emboss_path
      character*1000       :: promals_path
      character*1000       :: scwrl_path
      character*1000       :: loopy_path
      character*1000       :: loopy_copy ! 2022 v1.0
      character*1000       :: code_path
      character*66         :: string66
      character*66         :: strong66
      character*3          :: acrap1,acrap2,acrap3,acrap4,acrap5,acrap6
      character*3          :: acrap7,acrap8
      character*1 aa(20)
      character*100        :: my_name
      character*200        :: my_directory
      character*1000       :: input_line(100)
      data aa /
     &   'A','R','N','D','C','Q','E','G','H','I',
     &   'L','K','M','F','P','S','T','W','Y','V' /

C BLOSUM62 matrix data taken from:
C https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

      integer blosum(20,20)
      data blosum /
     &   4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
     &  -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
     &  -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
     &  -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
     &   0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
     &  -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
     &  -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
     &   0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
     &  -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
     &  -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
     &  -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
     &  -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
     &  -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
     &  -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
     &  -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
     &   1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
     &   0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
     &  -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
     &  -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
     &   0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4/

      blah(1:26)='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      blah(27:52)='abcdefghijklmnopqrstuvwxyz'
      blah(53:78)='!@#$%^&*()_+-={}[]|\:;<>,.'

      dipep(1:20)='ACDEFGHIKLMNPQRSTVWY'

      num_peps=0

      call init_random_seed()

C initialize the stored input lines to spaces

      do n=1,100
        do m=1,1000
          input_line(n)(m:m)=' '
        enddo
      enddo

      num_het_typs=0 ! number of HET residue types to ignore

      do i=1,10000
        write(sequenceU(i),'(i1)')mod(i,10)
      enddo

C first, make sure we have a pristine directory by deleting any scwrl or
C loopy files...

      call system('rm scwrl* loop* my_identity')
      call system('whoami > my_identity')
      open(unit=11,file='my_identity',status='unknown')
      read(11,'(a100)')my_name
      close(11)
      call system('pwd > my_directory')
      open(unit=11,file='my_directory',status='unknown')
      read(11,'(a200)')my_directory
      close(11)
      
C note that we assume the first argument is the cutoff score for
C determining a correct alignment between sequence and pdb

      ninput_lines=0
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000

C v4.6 also read "yes" or "no" for separate_hetatms
C v4.6 also read "yes" or "no" for use_norestrict
C v4.7 also read sidechain_mode as an integer

ccc   read(string1000,*)final_name ! final name for pdb
      read(string1000,*)final_name,char3,char3b,
     &                  sidechain_mode  ! =0 for all-atom 

      separate_hetatms=.false.
      if(char3.eq."yes") then
        separate_hetatms=.true.
        write(*,*)
        write(*,*)'will write fixed HETATMs to separate file'
        write(*,*)
      endif

      use_norestrict=.false.
      if(char3b.eq."yes") then
        use_norestrict=.true.
        write(*,*)
        write(*,*)'will use unrestricted confs for dipeptides'
        write(*,*)
      endif

      do lll=1,1000 
        sys_string9(lll:lll)=' '
        sys_string9b(lll:lll)=' '
        sys_string9c(lll:lll)=' '
      enddo
      sys_string9(1:len(trim(final_name)))=trim(final_name)
      sys_string9b(1:len(trim(final_name)))=trim(final_name)
      sys_string9c(1:len(trim(final_name)))=trim(final_name)
      k=len(trim(final_name))
      sys_string9(k+1:k+13)='.AHEMODEL.pdb'
      sys_string9b(k+1:k+12)='.CARTOON.pdb'
      sys_string9c(k+1:k+12)='.HETATMs.pdb'
      write(*,*)'final file 1 will be called ',sys_string9(1:k+13)
      write(*,*)'final file 2 will be called ',sys_string9b(1:k+12)
      if(separate_hetatms) then
        write(*,*)'final file 3 will be called ',sys_string9c(1:k+12)
      endif

C 2023 read the skip logicals

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,*)char3b,         ! skip_scwrl
     &                  char3c,         ! skip_dssp
     &                  char3d          ! skip_loopy

      skip_scwrl=.false.
      if(char3b.eq."yes") then
        skip_scwrl=.true.
        write(*,*)
        write(*,*)'will skip scwrl - assumes intra-template complete'
        write(*,*)
      endif

      skip_dssp=.false.
      if(char3c.eq."yes") then
        skip_dssp=.true.
        write(*,*)
        write(*,*)'will skip dssp'
        write(*,*)
      endif

      skip_loopy=.false.
      if(char3d.eq."yes") then
        skip_loopy=.true.
        write(*,*)
        write(*,*)'will skip loopy - assumes rigid domains are complete'
        write(*,*)
      endif

C read the path to the dipeptide pdb files

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')dipeptide_path
      dipeptide_path=adjustl(dipeptide_path)
      klm=len(trim(dipeptide_path))
ccc   write(*,*)'klm = ',klm

C v3.0 read the path to the ss2 files
C 2022 v1.0 ss2_path now defunct

c     read(5,'(a1000)')string1000
c     ninput_lines=ninput_lines+1
c     input_line(ninput_lines)=string1000
c     read(5,'(a1000)')string1000
c     ninput_lines=ninput_lines+1
c     input_line(ninput_lines)=string1000
c     read(5,'(a1000)')string1000
c     ninput_lines=ninput_lines+1
c     input_line(ninput_lines)=string1000
c     read(string1000,'(a1000)')ss2_path
c     ss2_path=adjustl(ss2_path)
c     klm2=len(trim(ss2_path))

C 2022 v1.0 read paths in this order:

C (1) dssp (2) emboss (3) promals (4) scwrl (5) loopy (6) localcode

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')dssp_path
      dssp_path=adjustl(dssp_path)
      len_dssp=len(trim(dssp_path))

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')emboss_path
      emboss_path=adjustl(emboss_path)
      len_emboss=len(trim(emboss_path))

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')promals_path
      promals_path=adjustl(promals_path)
      len_promals=len(trim(promals_path))

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')scwrl_path
      scwrl_path=adjustl(scwrl_path)
      len_scwrl=len(trim(scwrl_path))

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')loopy_path
      loopy_path=adjustl(loopy_path)
      len_loopy=len(trim(loopy_path))

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,'(a1000)')code_path
      code_path=adjustl(code_path)
      len_code=len(trim(code_path))

C now read a bunch of other parameters

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,*)score_cut,rms_crit,clash_dist1,clash_dist2,
     &                  nfails_chn_max,nfails_res_max,acrap1

C 2023_v1.1

      clash_max=max(clash_dist1,clash_dist2)

C v1.6 - allow us to continue if a chain in the template has no match

      if(acrap1(1:3).eq.'yes') then
        continue_regardless=.true.
        write(*,*)'continue_regardless = TRUE'
      else
        continue_regardless=.false.
        write(*,*)'continue_regardless = FALSE'
      endif    

      write(*,*)'using an alignment cutoff score of ',score_cut
      write(*,*)'using an RMS cutoff for dipeptides of ',rms_crit
      write(*,*)'using an ATOM-ATOM clash distance of ',clash_dist1
      write(*,*)'using an ATOM-HETATM clash distance of ',clash_dist2
      write(*,*)'# of tail failures before extending = ',nfails_chn_max
      clash_aa2=clash_dist1**2
      clash_ha2=clash_dist2**2

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000

      write(*,*)'HEADS-UP!!!'
      write(*,*)'if an error follows immediately then you may be '
      write(*,*)'using an old input file (pre v1.3) - you will '
      write(*,*)'need to add a final word "yes" or "no" to the '
      write(*,*)'end of the line beginning with "align_mode"'
      write(*,*)'you probably want a "no" to indicate no truncation'
      write(*,*)'of Nter or Cter tails'

      read(string1000,*)align_mode,acrap1,acrap2,acrap3,acrap4,acrap5,
     &                             acrap6 ! added for v3.6

C align_mode is used to determine behavior for needle alignments
C if = 1 then do alignments and proceed with building model
C    = 2 then do alignments and stop
C    = 3 then don't do alignments - read them from file

      i_go_up=0 ! default setting - not strictly necessary
      if(acrap1(1:3).eq.'yes') then
        membrane_dive=.true.
        write(*,*)'membrane_dive = TRUE'
      else
        membrane_dive=.false.
        write(*,*)'membrane_dive = FALSE'
      endif    

      if(acrap2(1:3).eq.'yes') then
        use_SS_constraints=.true.
        write(*,*)'use_SS_constraints = TRUE'
      else
        use_SS_constraints=.false.
        write(*,*)'use_SS_constraints = FALSE'
      endif    

      write(*,*)
      write(*,*)'possible dssp types (in priority order) are:'
      write(*,*)
      write(*,*)'H = alpha-helix'
      write(*,*)'B = beta-bridge'
      write(*,*)'E = beta-sheet'
      write(*,*)'G = 310-helix'
      write(*,*)'I = pi-helix'
      write(*,*)'T = H-bonded turn'
      write(*,*)'S = bend'
      write(*,*)'C = coil - not seen it yet, may not be used'
      write(*,*)

      num_dssp_types=len(trim(acrap3))
      do n=1,num_dssp_types
        my_dssp_type(n)=acrap3(n:n)
      enddo
      do n=1,num_dssp_types
        write(*,*)'we will look for dssp SS type = ',my_dssp_type(n)
      enddo

      if(acrap4(1:3).eq.'yes') then
        rebuild_tails=.true.
        rebuild_auto=.false.
        write(*,*)'rebuild_tails = TRUE'
      elseif(acrap4(1:3).eq.'no') then
        rebuild_tails=.false.
        rebuild_auto=.false.
        write(*,*)'rebuild_tails = FALSE'
      elseif(acrap4(1:3).eq.'aut') then
        rebuild_tails=.false.
        rebuild_auto=.true.
        write(*,*)'rebuild_tails = FALSE'
        write(*,*)'rebuild_auto  = TRUE'
      endif    

C there are three options for truncate_tails - "yes" means we need to read
C a file with the definitions in - "all" means build no N-terminal and no
C C-terminal tails...

      if(acrap5(1:3).eq.'yes') then
        truncate_tails=.true.
        build_no_tails=.false.
        write(*,*)'truncate_tails = TRUE'
      elseif(acrap5(1:3).eq.'no') then
        truncate_tails=.false.
        build_no_tails=.false.
        write(*,*)'truncate_tails = FALSE'
      elseif(acrap5(1:3).eq.'all') then
        truncate_tails=.false.
        build_no_tails=.true. 
        write(*,*)'build_no_tails = TRUE'
      endif    

      if(acrap6(1:3).eq.'yes') then
        skip_all_loops=.true.
        write(*,*)'skip_all_loops = TRUE - will not build any loops'
        write(*,*)
        write(*,*)'note that this should probably be combined with '
        write(*,*)'truncate_tails = all'
        write(*,*)
      else
        skip_all_loops=.false.
        write(*,*)'skip_all_loops = FALSE - will build all loops'
      endif

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,*)nmodels_desired,nloopy_confs,
     &                  min_loop_size,acrap4,acrap5,
     &                  nrelax,frelax ! v4.4

      if(acrap4(1:3).eq.'yes') then
        truncate_template=.true.
        write(*,*)'truncate_template = TRUE'
      else
        truncate_template=.false.
        write(*,*)'truncate_template = FALSE'
      endif    

      if(acrap5(1:3).eq.'yes') then
        rosetta_symmetry=.true.
        write(*,*)'rosetta_symmetry = TRUE'
      else
        rosetta_symmetry=.false.
        write(*,*)'rosetta_symmetry = FALSE'
      endif    

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000

C v4.1 add in iskiphetcheck to allow skipping of agglomeration step when
C selecting HETATM entries to retain
C put in a write statement to warn users that this is a late addition

      write(*,*)
      write(*,*)'if we crash next then add iskiphetcheck=0 or 1'
      write(*,*)'on the line beginning use_promals3D'
      write(*,*)'put "1" to include all HETATMs'
      write(*,*)'put "0" to iteratively include only close HETATM'
      write(*,*)
      write(*,*)'and add nbacktrack=0 or greater immediately after'
      write(*,*)'and nbacktries also...'
      write(*,*)'and iworkoutwards also'
      write(*,*)
      read(string1000,*)char4,acrap7,ss_weight,id_thr,
     &                  helix_cut,sheet_cut,iskiphetcheck,
     &                  nbacktrack,nbacktries,char3

C v4.3 nbacktrack is # of residues to backtrack
C      nbacktries is # of times to use backtracked position before
C                 potentially modifying it by one residue at a time

C v4.4 add "yes" for iworkoutwards - this logical controls what happens
C when we repeatedly fail to build a tail - if "no" (default) then we
C work inwards, removing structure from the inner-most part of the
C protein - if "yes" then we work outwards, removing structure from the
C outer-most part of the protein: this could be useful for nascent chain
C building when the nascent chain has structure outside of the tunnel

      iworkoutwards=.false.
      if(char3.eq."yes") then
        iworkoutwards=.true.
        write(*,*)
        write(*,*)'WARNING iworkoutwards invoked'
        write(*,*)'this is best used for nascent chain building'
        write(*,*)
      endif

C v4.1 - here we'll allow users to enter "yes" "no" "2019" to select
C whether, or which version of promals3D to use... note that we replaced
C acrap6 with char4 in the above read command

      use_promals3D=.false.
      use_promals3D_2019=.false.
      if(char4(1:3).eq.'yes') then
        use_promals3D=.true.
        write(*,*)'use_promals3D = TRUE'
        write(*,*)'we will use uniref90 that came with promals3D'
        write(*,*)
        write(*,*)'overriding use_SS_constraints & truncate_template'
        write(*,*)
        use_SS_constraints=.false.
        truncate_template=.false.
      elseif(char4(1:4).eq.'2019') then
        use_promals3D_2019=.true.
        write(*,*)'use_promals3D_2019 = TRUE'
        write(*,*)'we will use uniref90_2019 version'
        write(*,*)
        write(*,*)'overriding use_SS_constraints & truncate_template'
        write(*,*)
        use_SS_constraints=.false.
        truncate_template=.false.
      else
        write(*,*)'use_promals3D      = FALSE'
        write(*,*)'use_promals3D_2019 = FALSE'
      endif

      if(acrap7(1:3).eq.'yes') then
        use_tmalign=.true.
        write(*,*)
        write(*,*)'use_tmalign = TRUE'
        write(*,*)'using 3D information in promals3D'
        write(*,*)
      else
        use_tmalign=.false.
        write(*,*)
        write(*,*)'use_tmalign = FALSE'
        write(*,*)'ignoring 3D information in promals3D'
        write(*,*)
      endif

C v3.9 - assume that we have no ss2 files until we successfully read one

      have_ss2_file=.false.

C set score_weight (weight given to sequence) so that they add up to 1

      score_weight=1.0-ss_weight

C we limit ourselves to 10 confs only because I'm too lazy to rewrite
C the commands that look for loopy output to allow the possibility of
C one, two or three digit pdb numbers
C in practice my guess is that if loopy is having trouble finding a conf
C that works with 10 candidates it'll have trouble also with 20...

      if(nloopy_confs.gt.10) then
        write(*,*)
        write(*,*)'nloopy_confs must be <=10 currently'
        write(*,*)'quitting now :('
        write(*,*)
        stop
      endif

C read strong1000 to find the input .pdb file name - the actual pdb file
C name is figured out later

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(string1000,*)strong1000

C read string1000 to find all words - each will be a .fasta...

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000

      nl=len(trim(string1000))
      write(*,*)'nl = ',nl
      nwords=0
      if(string1000(1:1).eq.' ') then
        ingap=1
        inword=0
      else
        ingap=0
        inword=1
        ibeg_word(nwords)=1
      endif
      do n=2,nl
        if(string1000(n:n).eq.' ') then
          if(inword.eq.1) then
            ingap=1
            inword=0
            iend_word(nwords)=n-1
          endif
        elseif(string1000(n:n).ne.' ') then
          if(ingap.eq.1) then
            ingap=0
            inword=1
            nwords=nwords+1
            ibeg_word(nwords)=n
          endif
        endif
      enddo
      iend_word(nwords)=nl  

      write(*,*)'just read ',nwords,' fasta entries'
      num_seqs=nwords
      num_pdbs=1 ! hardwired to 1 - don't see any use for multiple pdbs

      write(*,*)'num_seqs = ',num_seqs
      write(*,*)'num_pdbs = ',num_pdbs

      write(*,*)
      write(*,*)'note that code is hardwired to allow a max '
      write(*,*)'of 500 chains per pdb file - seems like overkill'
      write(*,*)

      allocate(fname_seq(1:num_seqs))
      allocate(nname_seq(1:num_seqs)) ! # characters in name
      allocate(seq_of_fasta(1:num_seqs))
      allocate(num_res_in_seq(1:num_seqs))
      allocate(fname_pdb(1:num_pdbs))
      allocate(fname_pdb_chn(1:num_pdbs,1:100000))
      allocate(fasta_pdb_chn(1:num_pdbs,1:100000))
      allocate(falign(1:num_seqs,1:num_pdbs,1:100000))
      allocate(falign2(1:num_seqs,1:num_pdbs,1:100000))
      allocate(ialign(1:num_seqs,1:num_pdbs,1:100000))
      allocate(salign(1:num_seqs,1:num_pdbs,1:100000))
      allocate(salign1(1:num_seqs,1:num_pdbs,1:100000))
      allocate(salign2(1:num_seqs,1:num_pdbs,1:100000))
      allocate(salign3(1:num_seqs,1:num_pdbs,1:100000))
      ialign=0 ! default to saying no seqs match to any pdb chain
      salign=0.0 ! default score
      salign1=0.0 ! default score
      salign2=0.0 ! default score
      salign3=0.0 ! default score

C assign the input .pdb name 

      fname_pdb(1)=trim(strong1000)

C now assign the .fasta names from the stems provided in the input_file

      do n=1,nwords
        write(*,*)'n & name ',n,string1000(ibeg_word(n):iend_word(n))
        write(*,*)'n ibeg iend ',n,ibeg_word(n),iend_word(n)
        nl=iend_word(n)-ibeg_word(n)+1

C v3.9 store length of the name of the filename stem for .ss2 reading

        nname_seq(n)=nl

        fname_seq(n)(1:nl)=string1000(ibeg_word(n):iend_word(n))
        fname_seq(n)(nl+1:nl+6)='.fasta'
        do m=nl+7,80
          fname_seq(n)(m:m)=' '
        enddo 
        fname_seq(n)=trim(fname_seq(n))
      enddo

C check to see if there are duplicate entries - this should be an error

      do n=1,nwords
        do m=n+1,nwords
          if(fname_seq(n).eq.fname_seq(m)) then
            write(*,*)
            write(*,*)'error: duplicate fasta sequence ',fname_seq(n)
            write(*,*)'remove the second copy from the input_file'
            write(*,*)'quitting now :('
            write(*,*)
            stop
          endif
        enddo
      enddo

C now download the .fasta files using wget...
C UPDATE - we first touch the file and then wc it to see if it already
C exists and contains data - if so, then we don't do the download

C v3.9 note that we could do this much more elegantly using the fortran
C intrinsic "inquire" now - no need to change for now...

      touch(1:6)='touch '
      wcount(1:3)='wc '
      fasta_get(1:36)='wget http://www.uniprot.org/uniprot/'
      do n=1,nwords
        nl=len(fname_seq(n))
        touch(7:7+nl-1)=fname_seq(n)
        call system(touch(1:7+nl-1))
        call system('rm wc.out')
        wcount(4:4+nl-1)=fname_seq(n)
        wcount(4+nl:4+nl+8)=' > wc.out'
        call system(wcount(1:4+nl+8))
        open(unit=95,file='wc.out',status='unknown')
        read(95,*)icount
        close(95)
        call system('rm wc.out')
        if(icount.eq.0) then
          fasta_get(37:37+nl-1)=fname_seq(n)
          fasta_get(37+nl:40+nl)=' -O '
          fasta_get(41+nl:41+nl+nl-1)=fname_seq(n)
          call system(fasta_get(1:41+nl+nl-1),istatus)
          write(*,*)
          write(*,*)"wget'ing fasta file # ",n,fname_seq(n),istatus
          write(*,*)
        else
          write(*,*)
          write(*,*)"NOT wget'ing fasta file # ",n,fname_seq(n)
          write(*,*)
        endif
ccc     write(*,*)'actual wget command = ',fasta_get(1:41+nl+nl-1)
ccc     write(*,*)
      enddo

C read any three-letter residue codes for HETNAM entries to ignore

      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000
      read(5,'(a1000)')string1000
      ninput_lines=ninput_lines+1
      input_line(ninput_lines)=string1000

      nl=len(trim(string1000))
      if(nl.eq.0) goto 101

      nwords=0
      if(string1000(1:1).eq.' ') then
        ingap=1
        inword=0
      else
        ingap=0
        inword=1
        ibeg_word(nwords)=1
      endif
      do n=2,nl
        if(string1000(n:n).eq.' ') then
          if(inword.eq.1) then
            ingap=1
            inword=0
            iend_word(nwords)=n-1
          endif
        elseif(string1000(n:n).ne.' ') then
          if(ingap.eq.1) then
            ingap=0
            inword=1
            nwords=nwords+1
            ibeg_word(nwords)=n
          endif
        endif
      enddo
      iend_word(nwords)=nl  

101   write(*,*)'just read ',nwords,' HETATM residues to omit'

      num_het_atms=0
      num_het_typs=0
      do n=1,nwords
        write(*,*)'n & name ',n,string1000(ibeg_word(n):iend_word(n))
        nl=iend_word(n)-ibeg_word(n)+1
        if(nl.ne.3) then
          write(*,*)'skipping entry without 3 characters ',
     &              string1000(ibeg_word(n):iend_word(n))
          cycle
        endif
        num_het_typs=num_het_typs+1
        nam_het_typs(num_het_typs)=
     &    string1000(ibeg_word(n):iend_word(n))
        write(*,*)'new HETNAM res to omit: ',nam_het_typs(num_het_typs)
      enddo

      do i=1,100000
        seq_tot_res(i:i)=' '
      enddo

C at this point we can open up any forced_alignment file and see if we
C have a forced assignment for ialign...

C v3.1 for ease of coding we will require that every pdb and chn have an
C assignment here - this will be checked later - we will used
C "nforced>0" as our "logical" flag to denote forced align was used

      nforced=0
      open(unit=95,file="forced_alignment",status='unknown')
841   read(95,*,end=842)i,j,k ! i=pdb#,j=chn#,k=seq#
      ialign(k,i,j)=1
      nforced=nforced+1
      goto 841
842   close(95)
      write(*,*)'WARNING - will override #alignments ',nforced

      do i=1,num_pdbs
        do j=1,500
          do k=1,80
            fasta_pdb_chn(i,j)(k:k)=' '
            fname_pdb_chn(i,j)(k:k)=' '
          enddo
        enddo
      enddo
      allocate(num_chns(1:num_pdbs))
      allocate(seq_of_pdb_chn(1:num_pdbs,1:10000))
      allocate(len_of_pdb_chn(1:num_pdbs,1:10000))

C next array is used to override len_of_pdb_chn when 
C use_SS_constraints = "yes" as this can go crazy when the template is
C longer than the target...
C v1.4

      allocate(len_of_pdb_chn_fin(1:num_pdbs,1:10000))

C set default values for membrane positions

      zmem_neg1=0.000
      zmem_neg2=0.000
      zmem_pos1=0.000
      zmem_pos2=0.000

C now we can loop over all pdb files and separate them into chains...
C note that we keep only "ATOM" and "HET/MSE" entries - and note that 
C if we find an alternative conformation (i.e. not " " or "A" in col 17)
C then we skip it also... *and* to avoid trouble later on we'll set the
C label "A" to " " (i.e. remove evidence of alternative conformer)

C update, we now also keep any HET entry whose residue name is not
C listed on the command line - we default to omitting water (HOH/SOL)
C and we write DUM residues to DUM_original.pdb only (we omit them from
C any clash-checking as we do also with HOH/SOL)

C v3.9 do a system call to remove any pre-existing HETERO_original.pdb:
C this has been shown to cause problems in the past...

      call system("rm HETERO_original.pdb",istatus)

      nremarks=0

      open(unit=21,file='HETERO_original.pdb',status='unknown')
      open(unit=22,file='DUM_original.pdb',status='unknown')
      do nnn=1,num_pdbs

        nch=0
        nres=999999
        chain='?'

C v1.5 keep track of all HETATM residues found in the template - some or
C all of these may be omitted depending on the input file

        num_hets_found=0

        fname_len=len(trim(fname_pdb(nnn)))-4
        write(*,*)'fname_len = ',fname_len
        write(*,*)'opening ',fname_pdb(nnn)

        open(unit=10,file=fname_pdb(nnn),status="unknown")

C v3.2 do a quick read of fname_pdb(nnn) to get #hetatms
C and use this to allocate memory for xh,yh,zh ; xi,yi,zi

        nat_tmp=0
3201    read(10,'(a20)',end=3202)char20
        if(char20(1:4).eq.'HETA') nat_tmp=nat_tmp+1
        if(char20(1:4).eq.'ATOM') then
          if(char20(18:20).eq.' DA'.or.
     &       char20(18:20).eq.' DC'.or.
     &       char20(18:20).eq.' DG'.or.
     &       char20(18:20).eq.' DT'.or.
     &       char20(18:20).eq.'A  '.or.
     &       char20(18:20).eq.'C  '.or.
     &       char20(18:20).eq.'G  '.or.
     &       char20(18:20).eq.'U  '.or.
     &       char20(18:20).eq.' RA'.or.
     &       char20(18:20).eq.' RC'.or.
     &       char20(18:20).eq.' RG'.or.
     &       char20(18:20).eq.' RU') then
            nat_tmp=nat_tmp+1
          endif
        endif
        goto 3201
3202    rewind(10)

C set a flag that indicates we haven't yet checked each HETATM to
C determine whether it should be kept...
C we set ih=1 for a HETATM if it has been put on the keep-me list

        allocate(ih(1:nat_tmp))
        allocate(ih_original(1:nat_tmp))
        ih=0
        ih_original=0

        allocate(sh(1:nat_tmp)) ! v4.5
        allocate(xh(1:nat_tmp))
        allocate(yh(1:nat_tmp))
        allocate(zh(1:nat_tmp))
        allocate(si(1:nat_tmp)) ! v4.5
        allocate(xi(1:nat_tmp))
        allocate(yi(1:nat_tmp))
        allocate(zi(1:nat_tmp))

C v4.5 mh/mi will store, for each of the HETATMs stored in
C HETERO_original.pdb, a 1 if we still want it written out at the end -
C if it ends up being part of a moving tail then it'll already be in
C loopy_inp.pdb so we don't need it, hence set to 0

        allocate(mh(1:nat_tmp)) ! v4.5
        allocate(mi(1:nat_tmp)) ! v4.5

C pre-v3.2 code follows

5       read(10,21,end=31)strong

C store all REMARK_AHE entries for writing out at the end...

        if(strong(1:10).eq.'REMARK_AHE') then
          nremarks=nremarks+1
          strang(nremarks)=strong
        endif

C note that we used to quit when reading ENDMDL statements as I was
C under the impression that this was only used for NMR structures.
C However it turns out that this is often used in biological assembly
C files, with the authors also not renaming the chains, so we need to be
C able to handle this situation instead - so from now on, we will
C require that all users manually remove additional models in NMR
C files, and remove all ENDMDL entries before rerunning

        if(strong(1:6).eq.'ENDMDL') then
          write(*,*)
          write(*,*)'WARNING!!!'
          write(*,*)'this may be a NMR structure - will assume that '
          write(*,*)'ENDMDL line means end of this model and will now'
          write(*,*)'proceed, ignoring all subsequent models'
          write(*,*)
          goto 31
        endif

        if(strong(1:3).eq.'ATO') then

C first we need to check if this is RNA or DNA - if it is then we
C instead add these to the list of HETATM entries - we do this because
C at the moment we have no way of treating them except as steric
C blockers - eventually we might be able to sample them too...

          if(strong(18:20).eq.' DA'.or.
     &       strong(18:20).eq.' DC'.or.
     &       strong(18:20).eq.' DG'.or.
     &       strong(18:20).eq.' DT'.or.
     &       strong(18:20).eq.'A  '.or.
     &       strong(18:20).eq.'C  '.or.
     &       strong(18:20).eq.'G  '.or.
     &       strong(18:20).eq.'U  '.or.
     &       strong(18:20).eq.' RA'.or.
     &       strong(18:20).eq.' RC'.or.
     &       strong(18:20).eq.' RG'.or.
     &       strong(18:20).eq.' RU') then
            num_het_atms=num_het_atms+1
            read(strong,729)char30,x,y,z
            sh(num_het_atms)=strong(1:66)
            mh(num_het_atms)=num_het_atms ! v4.5 assume we want it 
            xh(num_het_atms)=x
            yh(num_het_atms)=y
            zh(num_het_atms)=z
            strong(1:6)='HETATM'
            write(21,21)strong
            goto 5
          endif

C otherwise, we can proceed with this 'ATOM' - do some minor tweaks to
C the names and then jump ahead to line 15

          if(strong(17:17).eq.' ') goto 15
          if(strong(17:17).eq.'A') then   
            strong(17:17)=' '
            goto 15
          elseif(strong(17:17).ne.'A') then   
            goto 5 ! skip it
          endif

C I think the following lines should never get executed since we should
C have renamed all MSE residues prior to running ahemodel...

ccc     elseif(strong(1:3).eq.'HET'.and.strong(18:20).eq.'MSE') then
ccc       strong(1:6)='ATOM  '
ccc       if(strong(13:15).eq.'SE ') strong(13:15)=' SD'
ccc       strong(18:20)='MET'
ccc       if(strong(17:17).eq.' ') goto 15
ccc       if(strong(17:17).eq.'A') then   
ccc         strong(17:17)=' '
ccc         goto 15
ccc       elseif(strong(17:17).ne.'A') then   
ccc         goto 5 ! skip it
ccc       endif

C now skip if atom is HET/SOL/DUM or HET and on list of omit-me restypes
C otherwise, write the HETATM to HETERO_original.pdb...

        elseif(strong(1:6).eq.'HETATM') then

C v1.5 keep track of any new HETATM residue types - we can write them
C out to the screen if align_mode=2

          inew=1
          do i3=1,num_hets_found
            if(strong(18:20).eq.nam_hets_found(i3)) then
              inew=0
              exit
            endif
          enddo
          if(inew.eq.1) then
            num_hets_found=num_hets_found+1
            nam_hets_found(num_hets_found)=strong(18:20)
          endif
          
          if(strong(18:20).eq.'SOL') goto 5
          if(strong(18:20).eq.'HOH') goto 5

C if membrane_dive then we need to check for OPM atoms - these will tell
C us where the membrane starts and stops... - once we have that
C information we can just write the atom to DUM_original.pdb and skip

          if(membrane_dive) then
            if(strong(13:20).eq.' N   DUM') then
              read(strong,'(46x,f8.3)')zmem_neg
            elseif(strong(13:20).eq.' O   DUM') then
              read(strong,'(46x,f8.3)')zmem_pos
            endif
            if(zmem_neg1.eq.0.000) zmem_neg1=zmem_neg
            if(zmem_pos1.eq.0.000) zmem_pos1=zmem_pos
            if(zmem_neg2.eq.0.000.and.
     &     abs(zmem_neg-zmem_neg1).gt.100.0) zmem_neg2=zmem_neg
            if(zmem_pos2.eq.0.000.and.
     &     abs(zmem_pos-zmem_pos1).gt.100.0) zmem_pos2=zmem_pos
              
          endif

C we write OPM DUM atoms to DUM_original.pdb (not HETERO_original.pdb)

          if(strong(13:20).eq.' N   DUM'.or.
     &       strong(13:20).eq.' O   DUM') then
            write(22,21)strong
          endif
          if(strong(18:20).eq.'DUM') goto 5 ! now totally skip DUM 

          do n=1,num_het_typs
            if(strong(18:20).eq.nam_het_typs(n)) then
              write(*,*)'skipping HETATM ',strong(1:30)
              goto 5
            endif
          enddo
          write(21,21)strong

C also add the coords to a list of het atoms so we can use them to
C do clash checking when we build tails...

          num_het_atms=num_het_atms+1
          read(strong,729)char30,x,y,z
729       format(a30,3f8.3)
          sh(num_het_atms)=strong(1:66) ! v4.5
          mh(num_het_atms)=num_het_atms ! v4.5 assume we want it 
          xh(num_het_atms)=x
          yh(num_het_atms)=y
          zh(num_het_atms)=z
        endif

        goto 5 ! if we got here then we're not an atom we want

15      continue ! come here to continue processing new ATOM

        read(strong,77)ntemp
77      format(22x,i4)

C ask if this is a new chain - don't use TER statements - just ask if
C the chainID is different or the current res# is less than last one

C v2.6: now we run into an issue with DSSP - it seems to choke on
C chainpdb files where all atom numbers are above 100,000 (PDH_TEST) -
C so we'll now have a chain-specific atom numbering...

        if(strong(22:22).ne.chain.or.ntemp.lt.nres) then

          natm_this_chain=0 ! v2.6
          if(nch.gt.1) close(11)
          write(*,*)'new chain found ',strong(22:22),ntemp,nch
          nch=nch+1
         fname_pdb_chn(nnn,nch)(1:fname_len)=fname_pdb(nnn)(1:fname_len)
         fname_pdb_chn(nnn,nch)(fname_len+1:fname_len+8)='.XXX.pdb'
         fasta_pdb_chn(nnn,nch)(1:fname_len)=fname_pdb(nnn)(1:fname_len)
         fasta_pdb_chn(nnn,nch)(fname_len+1:fname_len+10)='.XXX.fasta'
          write(ajunk,88)nch
88        format(i3)
          if(ajunk(1:1).eq.' ') ajunk(1:1)='0'
          if(ajunk(2:2).eq.' ') ajunk(2:2)='0'
          if(ajunk(3:3).eq.' ') ajunk(3:3)='0'
          fname_pdb_chn(nnn,nch)(fname_len+2:fname_len+4)=ajunk
          fasta_pdb_chn(nnn,nch)(fname_len+2:fname_len+4)=ajunk
          write(*,*)'check me ',fname_pdb_chn(nnn,nch)
          open(unit=11,file=fname_pdb_chn(nnn,nch),status="unknown")
        endif
        chain=strong(22:22)
        nres=ntemp
ccc     if(strong(1:3).eq.'ATO'.or.strong(1:3).eq.'HET') 
ccc  &    write(11,21)strong

C v2.6 here's the new write statement...

        if(strong(1:3).eq.'ATO'.or.strong(1:3).eq.'HET') then
          natm_this_chain=natm_this_chain+1
          write(11,'(a4,i7,a69)')strong(1:4),natm_this_chain,
     &                           strong(12:80)
        endif
        last_strong=strong
        goto 5
21      format(a80)
31      close(10)
        close(11)

        num_chns(nnn)=nch
        write(*,*)'num_chns in ',nnn,' = ',num_chns(nnn)

      enddo ! go back and analyze the next .pdb file

C v4.5 get max # of chains in any pdb so we can be more sensible about
C allocating memory for iam_in_pdb_chn etc...

      max_chn=0
      do nnn=1,num_pdbs
        max_chn=max(max_chn,num_chns(nnn))
      enddo
      write(*,*)
      write(*,*)'max_chn = ',max_chn
      write(*,*)

C write out info if membrane_dive

      if(membrane_dive) then
        zmem_ave1=(zmem_neg1+zmem_pos1)/2.0
        zmem_ave2=(zmem_neg2+zmem_pos2)/2.0
        write(*,*)
        write(*,*)'after reading template pdb zmem_neg1 = ',zmem_neg1 
        write(*,*)'after reading template pdb zmem_pos1 = ',zmem_pos1 
        write(*,*)'after reading template pdb zmem_ave1 = ',zmem_ave1 
        if(zmem_neg2.ne.0.000) then

C before going further we need to be careful to make sure that zmem(1)
C is at lower z values than zmem(2) - this is not guaranteed if an OM
C DUM atom is read before an IM DUM atom - this should fix this issue

          if(zmem_neg2.lt.zmem_neg1) then
            ztmp_neg1=zmem_neg1
            ztmp_pos1=zmem_pos1
            ztmp_ave1=zmem_ave1
            zmem_neg1=zmem_neg2
            zmem_pos1=zmem_pos2
            zmem_ave1=zmem_ave2
            zmem_neg2=ztmp_neg1
            zmem_pos2=ztmp_pos1
            zmem_ave2=ztmp_ave1
          endif
          write(*,*)'after reading template pdb zmem_neg2 = ',zmem_neg2
          write(*,*)'after reading template pdb zmem_pos2 = ',zmem_pos2
          write(*,*)'after reading template pdb zmem_ave2 = ',zmem_ave2
        endif
        write(*,*)'these are OPM coords used to guide tails out of mem'
        write(*,*)
      endif

      close(21) ! close HETERO_original.pdb

      write(*,*)
      write(*,*)'total number of HETATM entries = ',num_het_atms
      write(*,*)

C read the fasta files for the input sequences

      do iii=1,num_seqs
        num_res_in_seq(iii)=0
        open(unit=21,file=fname_seq(iii),status='unknown')
c       open(unit=21,file=fname_seq(iii)(1:nname_seq(iii)),
c    &       status='unknown')
        read(21,*) ! skip top line
26      read(21,27,end=28)seq_chunk
27      format(a60)

        do i=1,len(trim(seq_chunk)) ! line don't need to have 60 chars
          if(seq_chunk(i:i).ne.' ') then

C v4.1 here is where we will change any non-canonical residues;
C the non-canonical one-letter codes are explained here:
C   http://ddbj.nig.ac.jp/ddbj/code-e.html

            if(seq_chunk(i:i).eq.'B') then     ! asn/asp
              write(*,*)'WARNING: changing res# ',num_res_in_seq(iii)+1,
     &                  ' in seq# ',iii,' from B to N'
              seq_chunk(i:i)='N'
            elseif(seq_chunk(i:i).eq.'Z') then ! gln/glu
              write(*,*)'WARNING: changing res# ',num_res_in_seq(iii)+1,
     &                  ' in seq# ',iii,' from Z to Q'
              seq_chunk(i:i)='Q'
            elseif(seq_chunk(i:i).eq.'X') then ! any amino acid
              write(*,*)'WARNING: changing res# ',num_res_in_seq(iii)+1,
     &                  ' in seq# ',iii,' from X to A'
              seq_chunk(i:i)='A'
            elseif(seq_chunk(i:i).eq.'J') then ! ile/leu
              write(*,*)'WARNING: changing res# ',num_res_in_seq(iii)+1,
     &                  ' in seq# ',iii,' from J to L'
              seq_chunk(i:i)='L'
            elseif(seq_chunk(i:i).eq.'U') then ! selenocysteine
              write(*,*)'WARNING: changing res# ',num_res_in_seq(iii)+1,
     &                  ' in seq# ',iii,' from U to C'
              seq_chunk(i:i)='C'
            endif

            num_res_in_seq(iii)=num_res_in_seq(iii)+1
            seq_of_fasta(iii)(num_res_in_seq(iii):num_res_in_seq(iii))=
     &                        seq_chunk(i:i)
          else 
            goto 28
          endif 
        enddo

        goto 26
28      continue

        close(21) ! v4.1 - can't believe this wasn't already added!

      enddo

C check the fastas
 
      do iii=1,num_seqs
        write(my_fmt,'(a,i0,a)'),'(a',num_res_in_seq(iii),')'
        write(*,*)'Important Check ',iii,num_res_in_seq(iii)
        write(*,my_fmt)seq_of_fasta(iii)(1:num_res_in_seq(iii))
      enddo

C v3.4 identify which dipeptides we need to read in...

      i_need_dipep=0
      do iii=1,num_seqs
        do jjj=1,num_res_in_seq(iii)-1
          do kkk=1,20
            if(seq_of_fasta(iii)(jjj:jjj).eq.    dipep(kkk:kkk)) iw=kkk 
            if(seq_of_fasta(iii)(jjj+1:jjj+1).eq.dipep(kkk:kkk)) jw=kkk 
          enddo
          i_need_dipep(iw,jw)=1
        enddo
      enddo

C now we need to extract fasta files from each pdb file...
C we will align these with the true fasta soon enough

      write(*,*)
      write(*,*)'IMPORTANT - getting fastas from SEQRES entries'
      write(*,*)

C this is the point that we return to if use_SS_constraints="yes" - we
C first proceed as if it's "no" then update some things..
C v1.4

      ifrst_time_thru=0
      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)
          len_of_pdb_chn_fin(nnn,mmm)=0
        enddo
      enddo

9119  continue
      ifrst_time_thru=ifrst_time_thru+1

      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)
          len_of_pdb_chn(nnn,mmm)=0
        enddo
      enddo
      maxlen=0
      do nnn=1,num_pdbs
        mmm=0
        nres=999999
        chainID_last='?'
        fname_len=len(trim(fname_pdb(nnn)))-4
        open(unit=10,file=fname_pdb(nnn),status="unknown")
105     read(10,21,end=106)strong
        if(strong(1:6).ne.'SEQRES') goto 105
        chainID=strong(12:12)
        if(chainID.ne.chainID_last) then
          mmm=mmm+1
          chainID_last=chainID
        endif
        do mno=1,13
          jbeg=(mno-1)*4+20
          jend=(mno-1)*4+22
          res_n=strong(jbeg:jend)
          chain='?'
          if((res_n).eq.'   ') cycle
          if((res_n).eq.'ALA') chain='A'
          if((res_n).eq.'ARG') chain='R'
          if((res_n).eq.'ASN') chain='N'
          if((res_n).eq.'ASP') chain='D'
          if((res_n).eq.'ASH') chain='D'
          if((res_n).eq.'CYS') chain='C'
          if((res_n).eq.'CYM') chain='C'
          if((res_n).eq.'GLN') chain='Q'
          if((res_n).eq.'GLU') chain='E'
          if((res_n).eq.'GLH') chain='E'
          if((res_n).eq.'GLY') chain='G'
          if((res_n).eq.'HIS') chain='H'
          if((res_n).eq.'ILE') chain='I'
          if((res_n).eq.'LEU') chain='L'
          if((res_n).eq.'LYS') chain='K'
          if((res_n).eq.'LYN') chain='K'
          if((res_n).eq.'MLY') chain='K'
          if((res_n).eq.'MET') chain='M'
          if((res_n).eq.'MSE') chain='M'
          if((res_n).eq.'PHE') chain='F'
          if((res_n).eq.'PRO') chain='P'
          if((res_n).eq.'SER') chain='S'
          if((res_n).eq.'THR') chain='T'
          if((res_n).eq.'TRP') chain='W'
          if((res_n).eq.'TYR') chain='Y'
          if((res_n).eq.'VAL') chain='V'
          if((res_n).eq.'XXX') chain='X' ! v1.6
          if(chain.eq.'?') then
            write(*,*)
            write(*,*)'SHIT - found a weird residue in '
            write(*,*)'file name = ',fname_pdb_chn(nnn,mmm)
            write(*,*)'offending residue is ',res_n
            write(*,*)'quitting now :('
            write(*,*)
            stop
          endif

C here's where we stop reading the sequence if, on the third time
C through, with use_SS_constraints, we're past the *desired* end of the
C template chain (this may be much shorter than the *actual* end)
C note, however, that we only do this if truncate_template = true
C v1.4

          if(truncate_template) then
            if(use_SS_constraints.and.ifrst_time_thru.eq.3.and.
     &         len_of_pdb_chn(nnn,mmm)+1.gt.len_of_pdb_chn_fin(nnn,mmm))
     &      cycle
          endif

          len_of_pdb_chn(nnn,mmm)=len_of_pdb_chn(nnn,mmm)+1
          seq_of_pdb_chn(nnn,mmm)(len_of_pdb_chn(nnn,mmm):
     &                            len_of_pdb_chn(nnn,mmm))=chain
        enddo ! now look at next of the 13 entries on each line
        maxlen=max(len_of_pdb_chn(nnn,mmm),maxlen)
        goto 105 ! go back and read another line of the pdb file
106     close(10)
      enddo

C v2.0 we may need a few extra residues to be added to maxlen to account
C for possible differences between the fasta submitted to ahemodel and
C the fast submitted by Tim when he calculated .ss2...
C     maxlen=maxlen+50

C v2.6 - UPDATE - let's instead make maxlen be the biggest of maxlen as
C defined above and the length of the query fastas...

      do iii=1,num_seqs
        maxlen=max(maxlen,num_res_in_seq(iii))
      enddo

      write(*,*)
      write(*,*)'maxlen = ',maxlen
      write(*,*)

      if(use_SS_constraints.and.ifrst_time_thru.ge.2) then
        deallocate(iam_in_pdb_chn)
        deallocate(iam_in_SS_elem)
        deallocate(my_DSSP_design)
        deallocate(my_PSIPRED_prd)
        deallocate(my_PSIPRED_cnf)
        deallocate(my_PSIPRED_vl1)
        deallocate(my_PSIPRED_vl2)
        deallocate(my_PSIPRED_vl3)
      endif

C v4.5 4 April 2021 - made the allocation statements much more sensible
C here - note that the system calls can fail later if we run out of
C memory - but it won't tell us that...

ccc   allocate(iam_in_pdb_chn(1:num_pdbs,1:100000,1:maxlen),stat=ist)
      allocate(iam_in_pdb_chn(1:num_pdbs,1:max_chn,1:maxlen),stat=ist)
      write(*,*)'ISTAT001 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate iam_in_pdb_chn - quitting'
        write(*,*)
        stop
      endif
      allocate(iam_in_SS_elem(1:num_pdbs,1:max_chn,1:maxlen),stat=ist)
      write(*,*)'ISTAT002 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate iam_in_SS_elem - quitting'
        write(*,*)
        stop
      endif
      allocate(my_DSSP_design(1:num_pdbs,1:max_chn,1:maxlen),stat=ist)
      write(*,*)'ISTAT003 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_DSSP_design - quitting'
        write(*,*)
        stop
      endif
ccc   allocate(my_PSIPRED_prd(1:num_seqs,1:max_chn,1:maxlen),stat=ist)
      allocate(my_PSIPRED_prd(1:num_seqs,1:num_pdbs,1:maxlen),stat=ist)
      write(*,*)'ISTAT004 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_PSIPRED_prd - quitting'
        write(*,*)
        stop
      endif
      allocate(my_PSIPRED_cnf(1:num_seqs,1:num_pdbs,1:maxlen),stat=ist)
      write(*,*)'ISTAT005 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_PSIPRED_cnf - quitting'
        write(*,*)
        stop
      endif
      allocate(my_PSIPRED_vl1(1:num_seqs,1:num_pdbs,1:maxlen),stat=ist)
      write(*,*)'ISTAT006 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_PSIPRED_vl1 - quitting'
        write(*,*)
        stop
      endif
      allocate(my_PSIPRED_vl2(1:num_seqs,1:num_pdbs,1:maxlen),stat=ist)
      write(*,*)'ISTAT007 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_PSIPRED_vl2 - quitting'
        write(*,*)
        stop
      endif
      allocate(my_PSIPRED_vl3(1:num_seqs,1:num_pdbs,1:maxlen),stat=ist)
      write(*,*)'ISTAT008 ',ist
      if(ist.ne.0) then
        write(*,*)
        write(*,*)'failed to allocate my_PSIPRED_vl3 - quitting'
        write(*,*)
        stop
      endif

      iam_in_pdb_chn=-1  ! assume no coords found for residues in pdb
      iam_in_SS_elem=0   ! assume no SS elements found for residues in pdb
      my_DSSP_design=' ' ! assume no DSSP designation for residues in pdb
      my_PSIPRED_prd=' ' ! assume no PSIPRED prediction for residues in pdb
      my_PSIPRED_cnf=' ' ! assume PSIPRED confidence is zero for residues in pdb
      my_PSIPRED_vl1=0.0 ! assume PSIPRED confidence is zero for residues in pdb
      my_PSIPRED_vl2=0.0 ! assume PSIPRED confidence is zero for residues in pdb
      my_PSIPRED_vl3=0.0 ! assume PSIPRED confidence is zero for residues in pdb

C having finished with this pdb file we can write out fastas...

      write(*,*)'check maxlen = ',maxlen
ccc   write(*,*)'check len of = ',len_of_pdb_chn(1,1)

      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)
          open(unit=12,file=fasta_pdb_chn(nnn,mmm),status='unknown')
          write(12,'(a14)')'>protein.fasta'
          write(my_fmt,'(a,i0,a)'),'(a',len_of_pdb_chn(nnn,mmm),')'
          write(12,my_fmt)seq_of_pdb_chn(nnn,mmm)
     &                 (1:len_of_pdb_chn(nnn,mmm))
          close(12)
        enddo
      enddo

      write(*,*)
      write(*,*)'IMPORTANT - now determining which residues are'
      write(*,*)'actually present in the pdb... this will set '
      write(*,*)'iam_in_pdb_chn(nnn,mmm,nres)'
      write(*,*)
      write(*,*)'we will also check the res#s here to make sure '
      write(*,*)'they match the SEQRES entries'
      write(*,*)

      do nnn=1,num_pdbs

C v1.5 allow a few mismatches per chain between SEQRES and ATOM entries
C in the template - if there's a mismatch then correct the SEQRES...

        nmismatch=0
        do mmm=1,num_chns(nnn)
          ntmp=0
          open(unit=11,file=fname_pdb_chn(nnn,mmm),status='unknown')
43        read(11,21,end=42)strong
          if(strong(13:16).eq.' N  ') then

C for cases where we're on the second time through with
C use_SS_constraints we need to make sure we don't read past the desired
C end of the template - this should be safe to add for all cases

            ntmp=ntmp+1
            if(ntmp.gt.maxlen) cycle

            read(strong(23:26),'(i4)')nres
            read(strong(55:60),'(f6.2)')fres
            iam_in_pdb_chn(nnn,mmm,nres)=int(fres) ! stores domain #

            read(strong(18:20),'(a3)')res_n
            chain='?'
            if((res_n).eq.'   ') cycle
            if((res_n).eq.'ALA') chain='A'
            if((res_n).eq.'ARG') chain='R'
            if((res_n).eq.'ASN') chain='N'
            if((res_n).eq.'ASP') chain='D'
            if((res_n).eq.'ASH') chain='D'
            if((res_n).eq.'CYS') chain='C'
            if((res_n).eq.'CYM') chain='C'
            if((res_n).eq.'GLN') chain='Q'
            if((res_n).eq.'GLU') chain='E'
            if((res_n).eq.'GLH') chain='E'
            if((res_n).eq.'GLY') chain='G'
            if((res_n).eq.'HIS') chain='H'
            if((res_n).eq.'ILE') chain='I'
            if((res_n).eq.'LEU') chain='L'
            if((res_n).eq.'LYS') chain='K'
            if((res_n).eq.'LYN') chain='K'
            if((res_n).eq.'MET') chain='M'
            if((res_n).eq.'MSE') chain='M'
            if((res_n).eq.'PHE') chain='F'
            if((res_n).eq.'PRO') chain='P'
            if((res_n).eq.'SER') chain='S'
            if((res_n).eq.'THR') chain='T'
            if((res_n).eq.'TRP') chain='W'
            if((res_n).eq.'TYR') chain='Y'
            if((res_n).eq.'VAL') chain='V'
            if((res_n).eq.'XXX') chain='X' ! v1.6
            if(chain.ne.seq_of_pdb_chn(nnn,mmm)(nres:nres)) then
              nmismatch=nmismatch+1
              write(*,*)
              write(*,*)'WARNING!!!'
              write(*,*)'mismatch in file ',fname_pdb_chn(nnn,mmm)
              write(*,*)'for residue # ',nres,' in chain ',mmm
              write(*,*)'SEQRES entries say ',seq_of_pdb_chn(nnn,mmm)
     &                                       (nres:nres)
              write(*,*)'ATOM   entries say ',chain
              if(nmismatch.le.5) then
                seq_of_pdb_chn(nnn,mmm)(nres:nres)=chain
                write(*,*)'we will simply fix the SEQRES entry ' 
                write(*,*)'although this change is strictly unnecessary'
                write(*,*)'and continue for now...'
                write(*,*)
              else
                write(*,*)'we have now found 5 errors in this chain'
                write(*,*)'so you will need to fix the SEQRES entries'
                write(*,*)'before continuing - quitting now :('
                write(*,*)
                stop
              endif
            endif

          endif
          goto 43
42        close(11)
        enddo ! do another chain
      enddo   ! do another pdb

      write(*,*)'nres = ',nres

      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)
          do lll=1,len_of_pdb_chn(nnn,mmm)
            write(*,987)nnn,mmm,lll,seq_of_pdb_chn(nnn,mmm)(lll:lll),
     &                              iam_in_pdb_chn(nnn,mmm,lll)
987         format('GOOD PLACE TO STOP ',3i6,3x,a,3x,4i6)
          enddo
        enddo
      enddo

C before trying to do alignments we should call dssp for each chain

      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)

          do lll=1,1000
            sys_string1(lll:lll)=' '
          enddo

C 2022 v1.0 use dssp_path here
C UPDATE do this only if skip_dssp is false

          if(.not.skip_dssp) then

            sys_string1(1:len_dssp)=dssp_path(1:len_dssp)
            sys_string1(len_dssp+1:len_dssp+4)=' -i '
            l0=len_dssp+5
C           sys_string1(1:40)='/home/LAB/BIN/dssp-2.0.4-linux-amd64 -i '
            jjj=len(trim(fname_pdb_chn(nnn,mmm)))
            kkk=len(trim(fname_pdb_chn(nnn,mmm)))-3
            sys_string1(l0:l0+jjj-1)=fname_pdb_chn(nnn,mmm)(1:jjj)
            sys_string1(l0+jjj+0:l0+jjj+3)=' -o '
            sys_string1(l0+jjj+4:l0+jjj+4+kkk-1)=
     &                  fname_pdb_chn(nnn,mmm)(1:kkk)
            sys_string1(l0+jjj+4+kkk:l0+jjj+4+kkk+3)='dssp'

            write(*,*)'using system call here ',istatus
            call system(sys_string1(1:l0+jjj+4+kkk+3),istatus)

            if(nnn.eq.1.and.mmm.eq.1) write(*,*)'check call to dssp '
            write(*,*)sys_string1(1:l0+jjj+4+kkk+3),' stat = ',
     &                istatus

          endif

        enddo
      enddo

C now proceed through each sequence in turn and ask which pdb chain
C gives the best alignment for it

      do iii=1,num_seqs

        best_score=-999.999

        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)

C v3.1 here we can skip the alignment entirely if we've specified a
C forced_alignment and this iii,nnn,mmm combo is not assigned - this can
C help a lot when we have massive files (e.g. polysomes)

            if(nforced.gt.0.and.ialign(iii,nnn,mmm).eq.0) cycle

C reset my_res_dssp_type for all residues each time

            my_res_dssp_type=' ' ! I assume this does all elements

C 2022 v1.0 option to skip reading dssp file here

            if(skip_dssp) goto 601

C first we need to read the dssp file associated with this pdb file
C note that we do this way more times than necessary - we should be able
C to just do this once for each dssp file but that would mean storing
C stuff in memory - the cost associated with rereading the dssp file
C however is likely to be negligible so don't worry about this now

            do lll=1,1000
              sys_string1(lll:lll)=' '
              sys_string2(lll:lll)=' '
              sys_string3(lll:lll)=' '
              sys_string4(lll:lll)=' '
              sys_string5(lll:lll)=' '
              sys_string6(lll:lll)=' '
              sys_string7(lll:lll)=' '
              sys_string8(lll:lll)=' '
            enddo

            jjj=len(trim(fname_pdb_chn(nnn,mmm)))
            sys_string1(1:jjj-3)=fname_pdb_chn(nnn,mmm)(1:jjj-3)
            sys_string1(jjj-2:jjj+1)='dssp'
            kkk=len_trim(sys_string1)
cv4.5       write(*,*)'need to examine '
cv4.5       write(*,*)sys_string1(1:kkk)
            open(unit=11,file=sys_string1(1:kkk),status='unknown')

401         read(11,201)strong
201         format(a80)
            if(strong(1:12).ne.'  #  RESIDUE') goto 401

C if we got here then we're about to read the information for each amino
C acid - we want each amino acid type and secondary structure type

            last_ss_type='?'
            i_am_in_ss=0 ! use this to track if we're in SS or not
            num_ss=0

501         read(11,201,end=601)strong

C if there is an exclamation mark then there's a gap in the sequence
C here - we skip this line but first we reset last_ss_type to a crazy
C type so that we recognize the next SS type as a new element

            if(strong(14:14).eq.'!') then
              last_ss_type=';'
              i_am_in_ss=0
              goto 501
            endif
            read(strong,'(5x,i5)')nres

C now look at the secondary structure type - we accept this residue as
C as secondary structure element if its type matches with any of the
C num_dssp_types specified on the command line

C we later merge elements that have a gap of a single residue between
C them and, after merging, we ignore any elements of only 1 or 2 res

C note that we also store the DSSP assignment (regardless of whether
C it's a match) so that we can write it out in a sequence alignment for
C sanity checking purposes later - if we want to only store the DSSP
C assignment if it is a match to a desired type comment the following
C line and uncomment the one 5-6 lines below...

            curr_ss_type=strong(17:17)
            my_res_dssp_type(nres)=curr_ss_type

            i_do_this_one=0
            do n=1,num_dssp_types
              if(curr_ss_type.eq.my_dssp_type(n)) then
                i_do_this_one=1
ccc             my_res_dssp_type(nres)=curr_ss_type
              endif
            enddo
            if(i_do_this_one.eq.1) then
              if(i_am_in_ss.eq.0) then
                num_ss=num_ss+1
                i_am_in_ss=1
                ibeg_ss(num_ss)=nres
                iend_ss(num_ss)=nres
              else
                iend_ss(num_ss)=nres
              endif
            else
             i_am_in_ss=0
            endif
            last_ss_type=curr_ss_type

C now go back and read another amino acid

            goto 501

C come here when done reading the dssp file
C 2022 v1.0 also come here when skip_dssp is true

601         continue

C at this point I think we should merge SS elements that are spaced by
C only one amino acid - a helix with a bulge should probably be
C considered a single helix rather than two separate ones...

            do_me_still=1 ! assume we look at all SS elements
            num_ss_tmp=0
            do n=1,num_ss
              if(do_me_still(n).eq.0) cycle
              num_ss_tmp=num_ss_tmp+1
              ibeg_ss_tmp(num_ss_tmp)=ibeg_ss(n)
              iend_ss_tmp(num_ss_tmp)=iend_ss(n)
              m=1
801           continue

C look to see if the beginning of the next element is within 1AA of the
C end of the current element - if it is then we should amalgamate them

              if(ibeg_ss(n+m).eq.iend_ss_tmp(num_ss_tmp)+2) then
                iend_ss_tmp(num_ss_tmp)=iend_ss(n+m)
                do_me_still(n+m)=0 ! make sure we skip this one next time
                m=m+1 ! now look at the *next* element also
                goto 801
              endif
            enddo

C now copy these 'temporary' SS elements back into the regular arrays
C BUT note that we'll drop SS elements that contain only 1 or 2 AAs...

            num_ss=0
            do n=1,num_ss_tmp
              if(iend_ss_tmp(n)-ibeg_ss_tmp(n)+1.gt.2) then
                num_ss=num_ss+1
                ibeg_ss(num_ss)=ibeg_ss_tmp(n)
                iend_ss(num_ss)=iend_ss_tmp(n)
cv4.5           write(*,*)'SS element ',num_ss,ibeg_ss(num_ss),
cv4.5&                                         iend_ss(num_ss)
              endif
            enddo

C now store the SS# associated with each residue - we will write this
C to file imminently - note that we do not count the last amino acid of
C the SS element - this appears to be necessary to make the SS elements
C act properly as constraints during the needle alignments - until we
C have checked this properly we should label this TEMP TEMP TEMP

            nss=0
            do n=1,num_ss
ccc           do m=ibeg_ss(n),iend_ss(n)-1
              do m=ibeg_ss(n),iend_ss(n)
                nss(m)=n
              enddo
            enddo 

C now - if required - write these SS elements to SS_elements.txt so that
C they can be used as constraints during sequence alignment 

C note that we only do this on the second time through - i.e. *after* we
C have trimmed each sequence to its desired length...

            open(unit=99,file='SS_elements.txt',status='unknown')
            if(use_SS_constraints.and.ifrst_time_thru.ge.2) then
cv4.5         write(*,*)'writing SS constraints ',
cv4.5&                   len_of_pdb_chn(nnn,mmm)
              do lll=1,len_of_pdb_chn(nnn,mmm)
                write(99,'(2i8)')lll,nss(lll)
              enddo
            else
cv4.5         write(*,*)'NOT writing SS constraints'
cv4.5         write(99,*)
            endif
            close(99)

C having written them out for use in needle we can make sure that the
C SS elements are stored in iam_in_SS_elem - note that we make sure we
C record *all* residues that are in the SS element (including iend_ss)

            do n=1,num_ss
              do m=ibeg_ss(n),iend_ss(n)
                if(m.gt.maxlen) cycle ! this needed in new version
                iam_in_SS_elem(nnn,mmm,m)=n 
              enddo
            enddo 

C v1.4 store DSSP designation for all residues that match the desired
C designations (usually just H & E)

            do lll=1,len_of_pdb_chn(nnn,mmm)
              my_DSSP_design(nnn,mmm,lll)=my_res_dssp_type(lll)
            enddo

C now prepare to do the needle alignment - first make the very
C complicated string that issues the command - note that it assumes that
C needle is already in the PATH - AHE's modified source code is in:
C
C   /home/LAB/BIN/EMBOSS-6.6.0_AHE/needle
C
C and the executable has been copied to /usr/bin/local/needle

            sys_string1(1:len_emboss)=emboss_path(1:len_emboss)
            sys_string1(len_emboss+1:len_emboss+1)=' '
            l0=len_emboss+2

            jjj=l0+len(trim(fname_seq(iii)))-1
            kkk=l0+len(trim(fname_seq(iii)))+1+
     &             len(trim(fasta_pdb_chn(nnn,mmm)))

            sys_string1(l0:jjj)=fname_seq(iii)
            sys_string1(1+jjj:1+jjj)=' '
            sys_string1(2+jjj:kkk)=fasta_pdb_chn(nnn,mmm)
            sys_string1(kkk+1:kkk+28) ='-gapopen 10.0 -gapextend 0.5 '

            jjj=kkk+30
            kkk=jjj+len(trim(fname_seq(iii)))-6
            lll=len(trim(fname_seq(iii)))-6
            sys_string1(jjj:kkk)=fname_seq(iii)(1:lll)
            jjj_beg=jjj
            sys_string1(kkk:kkk)='_'
            jjj=kkk+1
            kkk=jjj+len(trim(fasta_pdb_chn(nnn,mmm)))-6
            lll=len(trim(fasta_pdb_chn(nnn,mmm)))-6
            sys_string1(jjj:kkk)=fasta_pdb_chn(nnn,mmm)(1:lll)
            sys_string1(kkk:kkk+4)='.out'
            jjj_end=kkk+4

C now read the alignment command to get the output filename in an easy
C to digest form - we'll call this falign...
C 2022 v1.0 the above way no longer works - I think because having
C slashes in the string causes read to go crazy - so now we store
C jjj_beg and jjj_end and use those to explicitly assign falign...

C           read(sys_string1(1:len(trim(sys_string1))),*)
C    &        a1,a2,a3,a4,a5,a6,a7,a8
C           falign(iii,nnn,mmm)=trim(a8)
            falign(iii,nnn,mmm)=sys_string1(jjj_beg:jjj_end)

c use the following for formatted write to file:
c           write(my_fmt,'(a,i0,a)'),'(a',len(trim(sys_string1)),')'
c           write(13,my_fmt)sys_string1(1:len(trim(sys_string1)))

C now make the call to needle if align_mode = 1 or 2 only...
C (if align_mode = 3 then we'll just read in the premade file - note
C  that this may well fail if we haven't previously run with = 1 or 2)

            if(align_mode.le.2)
     &        call system(sys_string1(1:len(trim(sys_string1))))
            if(align_mode.le.2) ! v3.5 might as well write this out too
     &        write(*,*)"system call: ",
     &                   sys_string1(1:len(trim(sys_string1)))

            sys_string2(1:13)='grep "Score" '
            sys_string2(14:14+len(falign(iii,nnn,mmm)))=
     &        falign(iii,nnn,mmm)
            kkk=14+len(falign(iii,nnn,mmm))+1
            sys_string2(kkk:kkk+9)='> temp.out'
            call system(sys_string2(1:kkk+9))
            open(unit=21,file='temp.out',status='unknown')
            read(21,*)a1,a2,f3
            close(21)

            salign(iii,nnn,mmm)=f3 ! store the alignment score
            if(truncate_template.and.ifrst_time_thru.eq.1) then
              salign1(iii,nnn,mmm)=f3 
            elseif(truncate_template.and.ifrst_time_thru.eq.2) then
              salign2(iii,nnn,mmm)=f3 
            elseif(truncate_template.and.ifrst_time_thru.eq.3) then
              salign3(iii,nnn,mmm)=f3 
            endif

          enddo
        enddo
      enddo

C before going on we want to look at every pdb/chain combo and decide
C which fasta it matches it best to - note that we'll only look at cases
C where the score exceeded the threshold - this should prevent us from
C building 'structures' for pdb/chains that don't have a good match...

C first, if we're using a forced alignment then we need to make sure
C that every pdb/chn got assigned to one and only one sequence

C v4.0 I no longer think that it should be fatal for a pdb/chn to have
C no assignment - instead we should continue - we'll catch these cases
C later if we have continue_regardless = no and stop then... so the
C lines below now represent a warning rather than an error...

      if(nforced.gt.0) then
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            iforced_assign=0
            do iii=1,num_seqs
              if(ialign(iii,nnn,mmm).eq.1) then
                iforced_assign=iforced_assign+1
              endif
            enddo
            if(iforced_assign.eq.0) then
              write(*,*)
              write(*,*)'you are using a forced_alignment file'
              write(*,*)'but pdb/chn ',nnn,mmm
              write(*,*)'has no forced assignment - this is a WARNING!'
              write(*,*)'to continue much further you will need to set '
              write(*,*)'continue_regardless = yes'
              write(*,*)
C             stop ! in v4.0 this is no longer a fatal error
            elseif(iforced_assign.gt.1) then
              write(*,*)
              write(*,*)'you are using a forced_alignment file'
              write(*,*)'but pdb/chn ',nnn,mmm
              write(*,*)'is assigned to multiple sequences'
              write(*,*)
              write(*,*)'this is a fatal error'
              write(*,*)
              write(*,*)'make sure that all pdb/chn combos are '
              write(*,*)'matched to *at most* one sequence'
              write(*,*)
              stop ! in v4.0 this remains a fatal error
            endif
          enddo
        enddo
      endif

      do nnn=1,num_pdbs
        do mmm=1,num_chns(nnn)

C check if this chain is already forced assigned to a sequence
C if it is then we can skip setting ialign for it, already done

          iforced_assign=0
          do iii=1,num_seqs
            if(ialign(iii,nnn,mmm).eq.1) then
              iforced_assign=1
              write(*,*)'chain/pdb ',nnn,mmm,
     &                  ' force-matches to seq ',iii
            endif
          enddo
          if(iforced_assign.eq.1) cycle

C for v1.6 we'll introduce best_score2 so that we can keep track of the
C best score obtained for each chain even if they don't exceed score_cut
C same with iii_best2...

          best_score=-999.999
          iii_best=-1
          best_score2=-999.999
          iii_best2=-1
          do iii=1,num_seqs
            if(salign(iii,nnn,mmm).gt.score_cut.and.
     &         salign(iii,nnn,mmm).gt.best_score) then
              best_score=salign(iii,nnn,mmm)
              iii_best=iii
            endif
            if(salign(iii,nnn,mmm).gt.best_score2) then
              best_score2=salign(iii,nnn,mmm)
              iii_best2=iii
            endif
          enddo
          if(iii_best.ne.-1) then
            write(*,*)'chain/pdb ',nnn,mmm,' matches to seq ',iii_best,
     &                ' with a score of ',best_score
            ialign(iii_best,nnn,mmm)=1
          else
            write(*,*)
            write(*,*)'WARNING!!!'
            write(*,*)'chain/pdb ',nnn,mmm,' does not match to a seq'
            write(*,*)'its best score was ',best_score2
            write(*,*)'with sequence # ',iii_best2
            write(*,*)
            if(continue_regardless) then
              write(*,*)'since continue_regardless = yes'
              write(*,*)'we will try to continue'
              write(*,*)'note that this means the atoms of this chain '
              write(*,*)'will be dropped from all further calculations'
              write(*,*)
              write(*,*)'if you want them to be present as steric'
              write(*,*)'occupiers of space you will probably need to'
              write(*,*)'make them HETATM entries instead - this will'
              write(*,*)'also require removing their SEQRES entries'
              write(*,*)'from the template_pdb'
              write(*,*)
            else
              write(*,*)'since continue_regardless = no'
              write(*,*)'this is fatal, so quitting :('
              write(*,*)
              stop
            endif
          endif
        enddo
      enddo

C v2.7 - I think we need to do a check now to make sure that the same
C sequence hasn't matched to chains of different types - there may be
C cases where this is desirable (can't think of any) but in general this
C should be a bad thing as it would suggest that we possibly have a
C heteroligomeric pdb but are only giving a single fasta...

C v4.0 - I think we should only do this if we are NOT using a
C forced_alignment - if we are, then this could undo the alignment...

      if(nforced.eq.0) then
        do iii=1,num_seqs

C first, go through and find the best score found for this sequence

          best_score=0.0
          nnn_best=0
          mmm_best=0

          do nnn=1,num_pdbs
            do mmm=1,num_chns(nnn)
              if(salign(iii,nnn,mmm).gt.best_score) then
                best_score=salign(iii,nnn,mmm)
                nnn_best=nnn
                mmm_best=mmm
              endif
            enddo
          enddo

C now reset ialign to zero for all matches that had lower scores

          do nnn=1,num_pdbs
            do mmm=1,num_chns(nnn)
              if(ialign(iii,nnn,mmm).eq.1.and.
     &           salign(iii,nnn,mmm).lt.best_score) then
                write(*,*)
                write(*,*)'WARNING!!!'
                write(*,*)'overriding match of chain/pdb ',nnn,mmm,
     &                    ' to seq# ',iii,' since this seq matches',
     &                    ' better to chain/pdb ',nnn_best,mmm_best
                write(*,*)'this usually indicates that the pdb is',
     &                    ' a heteroligomer and you are missing a fasta'
                write(*,*)
                ialign(iii,nnn,mmm)=0
              endif
            enddo
          enddo
        enddo
      endif
       
C we should also do a check to make sure that every input sequence found a
C match - if one didn't then we may have the wrong .fasta...

C v4.0 - this should definitely be done, regardless of forced_alignment
C - and it definitely should be a fatal error...

      do iii=1,num_seqs
        imatched=0
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            if(ialign(iii,nnn,mmm).eq.1) imatched=1
          enddo
        enddo
        if(imatched.eq.0) then
          write(*,*)
          write(*,*)'WARNING!!!'
          write(*,*)'sequence ',iii,' does not match to a pdb/chain'
          write(*,*)'you may have the wrong sequence or you may need'
          write(*,*)'to force an alignment for this sequence to a chain'
          write(*,*)'i.e. add the file forced_alignment'
          write(*,*)
          write(*,*)'alternatively, set score_cut to a lower value'
          write(*,*)'in the input file for ahemodel'
          write(*,*)
          write(*,*)'quitting :('
          write(*,*)
          stop
        endif
      enddo

C v3.1 we will also write out an alignment_summary file - this could be
C used as a forced_alignment file in subsequent submissions...

C v4.0 - write this fact to the screen so that we don't forget it!

      write(*,*)
      write(*,*)'writing summary of the final winning alignments to '
      write(*,*)'alignment_summary - this could be used as a '
      write(*,*)'forced_alignment file in subsequent submissions to '
      write(*,*)'ahemodel - this might save time with very large '
      write(*,*)'systems such as polysomes...'
      write(*,*)

      open(unit=95,file="alignment_summary",status='unknown')
      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            if(ialign(iii,nnn,mmm).eq.1) write(95,'(3i8)')nnn,mmm,iii
          enddo
        enddo
      enddo
      close(95)

C at this point we know, for each sequence, which pdb/chain it matches
C to - we can now go through and find out where the gaps and tails are

C first, we should figure out, for each sequence, which residues of
C which pdb/chain are exact matches, which are partial matches, and
C which are gaps/extensions - our goal should probably be to write out a
C structure for each pdb/chain for which there is a match - then
C amalgamate them into a single file, then run loopy and add tails...

      num_frags=0
      nres_scwrl=0

      num_tot_res=0 ! total number of residues in final model

C zero out the following arrays

      my_scwrl_res=0
      my_tot_res  =0
      i_done_res  =0 ! used for passing res's to loopy

      num_tot_loops=0
      num_tot_nters=0
      num_tot_cters=0
      num_act_cters=0 ! actual already-done C-terminu
      res_tot_res  =0 ! intra-chain res # to be written out

C zero out the net charge on the scwrl model

      qqq_tot_scwrl=0.0

C set i_am_original and i_old_unstruc to zero..
C i_am_original stores whether the residue was found in the pdb file -
C it doesn't care how good the match was, just whether there was a match
C i_old_unstruc stores whether the residue was model-built in the pdb
C file - it looks for beta=0 and only does this for identical matches,
C since non-identical matches will be model-built here anyway

      my_domain=0
      my_bfactr=-1.0 ! default value
      i_am_original=0
      i_old_unstruc=0
      i_am_by_loopy=0 ! assume not placed by loopy

      i_in_pdb_fasta_tot=0 ! set the global version to zero 
      i_in_dom_fasta_tot=-1 ! set the global version to -1 v2.0!!!
      i_in_SS_elemnt_tot=0 ! set the global version to zero 

      open(unit=23,file='scwrl_inp.pdb',status='unknown')
      open(unit=24,file='scwrl.fasta',status='unknown')
      open(unit=25,file='final_sequence_alignments.txt',
     &     status='unknown')

C v1.4 keep track of total # of chains for writing purposes and use this
C to index #of seq identities, positives, gaps for each alignment

      num_chns_tot=0
      chain_nam='XXXXXX' ! uniprot ID of each chain in final model

C remember that the following are all arrays dimensioned by #chains so
C it's appropriate to be zeroing them all outside of the main loop...

      nres_totl=0 ! #res in each chain (sanity check)
      nres_iden=0 ! #res in each chain that match exactly with template
      nres_posi=0 ! #res that match well but not exactly
      nres_alin=0 ! #res that align but don't match well
      nres_gaps=0 ! #res that have no match in template
      nres_befr=0 ! #res *before* alignment starts
      nres_aftr=0 ! #res *after* alignment starts
      nres_tmplte_gaps=0 ! #res in template that have no match in target
      nres_tmplte_befr=0 ! #res in template *before* alignment starts
      nres_tmplte_aftr=0 ! #res in template *after* alignment starts

C v2.1 the following used to assign code-local chainIDs to all fragments
C - these will not be used in the final write-out of final_pdb...

      num_chainID=1
      nres_this_chainID=0

C v2.9 keep counter of total res# for occ_psipred

      nres_final=0

C v2.7 open up a file storing the DSSP and PSIPRED correspondences

      open(unit=93,file='DSSP_vs_PSIPRED.txt',status='unknown')

      do iii=1,num_seqs

C v1.7 keep track of whether we already aligned this sequence - we can
C skip doing the (costly with promals3D) alignment again...

        ialreadydidalignment=0
        ididthisseq=0

        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)

            do n=1,100000
              sequenceA(n:n)=' '
              sequenceB(n:n)=' '
              sequenceM(n:n)=' '
              sequenceS(n:n)=' '
              sequenceD(n:n)=' '
              sequenceP(n:n)=' '
              sequenceC(n:n)=' '
              sequenceT(n:n)=' '
              sequenceU(n:n)=' '
            enddo

            i_am_in_a_del=0  ! assume no residues are in deletions

C note that i_in_pdb_fasta records whether a residue from the target
C matches to one for which we have coordinates in the template

C i_in_pdb_fasta=1 for any kind of match to another residue
C i_in_pdb_fasta=2 for similar residues
C i_in_pdb_fasta=3 for identical residues

C note that i_in_pdb_fasta is effectively dimensioned 1:nres1 where
C nres1 is the number of residues in the TARGET...

C i_in_SS_elemnt=1 means this res is in a desired SS element (H,E,etc)

            i_in_pdb_fasta=0 ! reset seq for each pdb  - this needs to 
                           ! change if/when we allow multiple templates
            i_in_SS_elemnt=0 ! reset seq for each pdb  - this needs to 
            
C v4.0 here, for sure, is the point at which we skip any atoms that
C belong to any unaligned sequences in the template pdb...
               
            if(ialign(iii,nnn,mmm).eq.0) cycle

C otherwise, we carry on and increment the seq and chain numbers done

            ididthisseq=ididthisseq+1
            num_chns_tot=num_chns_tot+1 ! v1.4
            chain_nam(num_chns_tot)=falign(iii,nnn,mmm)(1:6)

            i_found_frst_struct_res=0
            nblanks=0
            nres1=0
            nres2=0
            nres2_tru=0 ! residues that are truly present in template
            my_nres2=0 ! BUGFIX in v1.5 this was not being reset 
            my_nres2_tru=0 ! this is new to v1.5
            ntot=0 ! records all positions: residues/inserts/deletes
            num_frags=num_frags+1 ! needs change if we have mult templs
            num_res_in_frag(num_frags)=num_res_in_seq(iii)

C v2.1 figure out what the chainID should be for this particular
C fragment - we're going to reuse each chainID over and over until its
C residue count exceeds 9999 then move on to the next chainID...
C to achieve this, we'll introduce:

C num_chainID - the number of the current chainID (1=A,2=B,etc)
C nres_this_chainID - #res so far allocated to the current chainID
C cur_chainID - the current actual chainID of this frag (A,B,etc)
C fin_chainID - the final writeout chainID of this frag (A,B,etc)

C note that there are many cases below where chainIDs are set - we are
C going to decide for each of these cases whether they should be
C assigned the code-local chainID (cur_chainID) or whether they should
C be the final-output chainID (fin_chainID) - e.g.

C  chn_tot_res(num_tot_res:num_tot_res) - fin_chainID - DONE
C ichn_tot_nters(num_tot_nters)         - cur_chainID - DONE
C ichn_tot_cters(num_tot_cters)         - cur_chainID - DONE
C ichn_tot_loops(num_tot_loops)         - cur_chainID - DONE

C *****************************************************************
C to do: (1) finish other cases where chainID is set
C        (2) decide on a new residue numbering system that will work
C            for all cases - I think we should just need an "offset"
C            for each fragment/chain that converts between chain-local
C            numbering and the cur_chainID-local numbering system
C            - maybe it would be helpful to again use the names:
C              fin_offset & cur_offset? (fin_offset would be zero...)

C first, figure out what the final chainID will be for this fragment...

            ijk=mod(num_frags,chainIDmax)
            if(ijk.eq.0) ijk=chainIDmax
            fin_chainID(num_frags)=blah(ijk:ijk)
            fin_offset(num_frags)=0 ! alwys 0 - inclddd 2B complete

C now figure out if we need to increment the code-local chainID due to
C the residue count exceeding 9999...

            if(nres_this_chainID+num_res_in_frag(num_frags).gt.9999)
     &      then
              num_chainID=num_chainID+1
              if(num_chainID.gt.chainIDmax) then
                write(*,*)
                write(*,*)'sorry, now on num_chainID= ',num_chainID
                write(*,*)'this exceeds the compiled max of ',chainIDmax
                write(*,*)'you must have a *lot* of residues...'
                write(*,*)'quitting now :('
                write(*,*)
                stop
              endif
              nres_this_chainID=num_res_in_frag(num_frags)
            else
              nres_this_chainID=nres_this_chainID+
     &                          num_res_in_frag(num_frags)
            endif

C having figured out whether num_chainID has incremented we can set
C cur_chainID for this fragment - note that there is no need for a mod()
C statement her as num_chainID can never exceed chainIDmax (see above)

C note how we will use cur_offset moving forward: we take the
C chain-local residue number and add cur_offset to get the residue
C number as it will appear in all intermediate pdb files...

            cur_chainID(num_frags)=blah(num_chainID:num_chainID)
            cur_offset(num_frags)=nres_this_chainID-
     &                            num_res_in_frag(num_frags)

            write(*,*)'CHECK setting ',num_frags,
     &                                 num_res_in_frag(num_frags),
     &                                 nres_this_chainID,
     &                                 cur_chainID(num_frags),
     &                                 cur_offset(num_frags)

C v1.9 store, for each fragment (i.e. chain), the sequence that it
C matches to - we'll use this later for writing out loop seqs...

            my_frags_seq_num(num_frags)=iii

            if(num_frags.eq.1) then
              my_frst_res(1)=1
            else           
              my_frst_res(num_frags)=my_last_res(num_frags-1)+1
            endif
            my_last_res(num_frags)=my_frst_res(num_frags)+
     &                             num_res_in_frag(num_frags)-1
            my_totl_res(num_frags)=my_last_res(num_frags)-
     &                             my_frst_res(num_frags)+1

C v1.7 we could call promals3D here, then write separate code to convert
C the output alignment into needle format, then retain all of the
C following lines - that would be very clean...

C v4.1 allow use_promals3D or use_promals3D_2019...

            if(use_promals3D.or.use_promals3D_2019) then

              if(ialreadydidalignment.eq.1) then 
                write(*,*)'for sequence # ',iii
                write(*,*)'promals3D alignment already done so skip now'
                do lll=1,1000
                  pro_string0(lll:lll)=' '
                enddo
                ll2=len(trim(falign_tmp))
                pro_string0(1:3)='cp '
                pro_string0(4:4+ll2)=falign_tmp(1:ll2)
                jjj=5+ll2
                kkk=jjj+len(trim(fname_seq(iii)))-6
                lll=len(trim(fname_seq(iii)))-6
                pro_string0(jjj:kkk)=fname_seq(iii)(1:lll)
                pro_string0(kkk:kkk)='_'
                jjj=kkk+1
                kkk=jjj+len(trim(fasta_pdb_chn(nnn,mmm)))-6
                lll=len(trim(fasta_pdb_chn(nnn,mmm)))-6
                pro_string0(jjj:kkk)=fasta_pdb_chn(nnn,mmm)(1:lll)
                pro_string0(kkk:kkk+13)='.out.promals3D'
                write(*,*)
                write(*,*)"system call: ",pro_string0(1:kkk+13)
                call system(pro_string0(1:kkk+13))
                falign(iii,nnn,mmm)=pro_string0(5+ll2:kkk+13)
                goto 499
              endif

              jjj=4+len(trim(fname_seq(iii)))
              kkk=5+len(trim(fname_seq(iii)))+
     &              len(trim(fasta_pdb_chn(nnn,mmm)))

              do lll=1,1000
                pro_string1(lll:lll)=' '
                pro_string2(lll:lll)=' '
                pro_string3(lll:lll)=' '
                pro_string4(lll:lll)=' '
              enddo

C v3.6 following lines work well only if the fastas contain proper 
C carriage returns at the end, so use the subsequent lines instead

c             pro_string1(1:4)='cat '
c             pro_string1(5:jjj)=fname_seq(iii)
c             pro_string1(1+jjj:1+jjj)=' '
c             pro_string1(2+jjj:kkk)=fasta_pdb_chn(nnn,mmm)
c             pro_string1(kkk+1:kkk+16) =' > promals3D.inp'
c             write(*,*)"system call: ",pro_string1(1:kkk+16)
c             call system(pro_string1(1:kkk+16))

              pro_string1(1:5)='(cat '
              pro_string1(6:jjj+1)=fname_seq(iii)
              pro_string1(2+jjj:18+jjj)=' ; echo "" ; cat '
              pro_string1(19+jjj:17+kkk)=fasta_pdb_chn(nnn,mmm)
              pro_string1(kkk+18:kkk+34) =') > promals3D.inp'
              write(*,*)"system call: ",pro_string1(1:kkk+34)
              call system(pro_string1(1:kkk+34))

C 2022 v1.0 - code used to split here to read different databases
C associated with different versions of promals3D - this is irrelevant
C now that promals_path is explicitly defined in the input file

              pro_string2(1:len_promals)=promals_path(1:len_promals)
              l0=len_promals+1
              pro_string2(l0:l0+23)=' promals3D.inp -fast 0 '
              pro_string2(l0+24:l0+54)='-outfile promals3D.out -id_thr '
              write(char4,'(f4.2)')id_thr
              pro_string2(l0+55:l0+58)=char4
              pro_string2(l0+59:l0+70)=' -ss_weight '
              write(char4,'(f4.2)')ss_weight
              pro_string2(l0+71:l0+74)=char4
              pro_string2(l0+75:l0+89)=' -score_weight '
              write(char4,'(f4.2)')score_weight
              pro_string2(l0+90:l0+93)=char4

c             pro_string2(35:64)='promals promals3D.inp -fast 0 '
c             pro_string2(65:95)='-outfile promals3D.out -id_thr '
c             write(char4,'(f4.2)')id_thr
c             pro_string2( 96: 99)=char4
c             pro_string2(100:111)=' -ss_weight '
c             write(char4,'(f4.2)')ss_weight
c             pro_string2(112:115)=char4
c             pro_string2(116:130)=' -score_weight '
c             write(char4,'(f4.2)')score_weight
c             pro_string2(131:134)=char4

C if we are not using tmalign (i.e. if we are ignoring 3D information
C when running promals3D) then we need to tack on to the system call

              if(use_tmalign) then
                write(*,*)"system call: ",pro_string2(1:l0+93)
                call system(pro_string2(1:l0+93))
              else
                pro_string2(l0+94:l0+103)='-tmalign 0'
                write(*,*)"system call: ",pro_string2(1:l0+103)
                call system(pro_string2(1:l0+103))
              endif

              pro_string3(1:len_code)=code_path(1:len_code)
              l0=len_code+1
              pro_string3(l0:l0+26)='/convert_promals3D_out.exe '
              pro_string3(l0+27:l0+40)='promals3D.out '
              jjj=l0+41
              kkk=jjj+len(trim(fname_seq(iii)))-6
              lll=len(trim(fname_seq(iii)))-6
              pro_string3(jjj:kkk)=fname_seq(iii)(1:lll)
              pro_string3(kkk:kkk)='_'
              jjj=kkk+1
              kkk=jjj+len(trim(fasta_pdb_chn(nnn,mmm)))-6
              lll=len(trim(fasta_pdb_chn(nnn,mmm)))-6
              pro_string3(jjj:kkk)=fasta_pdb_chn(nnn,mmm)(1:lll)
              pro_string3(kkk:kkk+13)='.out.promals3D'
              write(*,*)"system call: ",pro_string3(1:kkk+13)
              call system(pro_string3(1:kkk+13))

C now update the name of the alignment file so that it ends with
C 'promals3D'...
C 2022 version - change 69 to l0+41 below so we don't have a space in
C filename that shouldn't be there...
 
              falign(iii,nnn,mmm)=pro_string3(l0+41:kkk+13)

C promals3D seems to throw an error when finding old files often
C so move the temporary promals3D blast folder to a different name

              pro_string4(1:23)='mv promals3D.inp_blast '
              jjj=24
              kkk=jjj+len(trim(fname_seq(iii)))-6
              lll=len(trim(fname_seq(iii)))-6
              pro_string4(jjj:kkk)=fname_seq(iii)(1:lll)
              pro_string4(kkk:kkk)='_'
              jjj=kkk+1
              kkk=jjj+len(trim(fasta_pdb_chn(nnn,mmm)))-6
              lll=len(trim(fasta_pdb_chn(nnn,mmm)))-6
              pro_string4(jjj:kkk)=fasta_pdb_chn(nnn,mmm)(1:lll)
              pro_string4(kkk:kkk+19)='.promals3D.inp_blast'
              write(*,*)"system call: ",pro_string4(1:kkk+19)
              call system(pro_string4(1:kkk+19))

              if(ialreadydidalignment.eq.0) then 
                falign_tmp=trim(adjustl(falign(iii,nnn,mmm)))
                ialreadydidalignment=1
              endif

C close the if for use_promals3D

            endif

499         write(*,*)'Now examining alignment: ',falign(iii,nnn,mmm)

            open(unit=11,file=falign(iii,nnn,mmm),status='unknown')
            indel=0 ! assume we're not in a deletion already

            write(*,*)'AHE HERE 001'
109         read(11,209,end=309)streng
209         format(a90)
            write(*,*)'AHE HERE 002',streng

            if(streng(1:1).eq.'#') goto 109

C put in a simple fix to correct problem with one line not having a '#'
C at the beginning of it - this means we skip it but not others like it

            nblanks=nblanks+1
            if(nblanks.eq.1) goto 109

C if we get here then we're starting to read the blocks of alignments

409         read(11,209)streng1 ! sequence of desired fasta
            if(streng1(1:1).eq.' ') goto 309 ! done
            read(11,209)streng0 ! identities & similarities
            read(11,209)streng2 ! sequence of template pdb

            write(*,*)streng1
            write(*,*)streng0
            write(*,*)streng2

            do i=1,50

C if we find a blank entry for the target we're at the end so stop

              if(streng1(21+i:21+i).eq.' ') goto 309

C otherwise we update the numbers of residues in target and template

              if(streng1(21+i:21+i).ne.'-') then
                nres1=nres1+1
                nres_totl(num_chns_tot)=nres_totl(num_chns_tot)+1 ! v1.4
              endif
              if(streng2(21+i:21+i).ne.'-') then
                nres2=nres2+1
              endif

C the next thing we do is find cases where a template residue that was
C used in the alignment is not actually present in the pdb - we set the
C alignment to nothing and remove its entry from the template sequence:
C v2.6 - can only do this if nres2>0 by definition....
      
              if(nres2.gt.0) then
                if(iam_in_pdb_chn(nnn,mmm,nres2).eq.-1) then
                  streng0(21+i:21+i)=' '
                  streng2(21+i:21+i)=' '
                endif
              endif

C v1.5 now we can increment the "true" number of res in the template - 
C we use this to determine template overhangs for stats purposes only

              if(streng2(21+i:21+i).ne.'-'.and.
     &           streng2(21+i:21+i).ne.' ') then
                nres2_tru=nres2_tru+1
              endif

C now we can go on and assess identities, gaps etc

              if(streng1(21+i:21+i).ne.'-') then
                if(streng0(21+i:21+i).eq.' ') then
                  nres_gaps(num_chns_tot)=nres_gaps(num_chns_tot)+1 ! v1.4
                endif
              else
                if(streng2(21+i:21+i).ne.'-'.and.
     &             streng2(21+i:21+i).ne.' ') then
                  nres_tmplte_gaps(num_chns_tot)=
     &            nres_tmplte_gaps(num_chns_tot)+1 ! v1.4
                endif
              endif
              if(streng2(21+i:21+i).ne.'-'.and.
     &           streng2(21+i:21+i).ne.' ') then
ccc           if(streng2(21+i:21+i).ne.'-') then
                if(iam_in_pdb_chn(nnn,mmm,nres2).ge.0.and.
     &             i_found_frst_struct_res.eq.0) then
                   nres_befr(num_chns_tot)=nres1-1
                   i_found_frst_struct_res=1
                endif
              endif
              if(streng0(21+i:21+i).eq.'|') then
                if(iam_in_pdb_chn(nnn,mmm,nres2).ge.0) then
                  nmatch=nmatch+1
                  i_in_pdb_fasta(nres1)=3 ! perfect match
                  i_in_pdb_fasta_tot(nres1+my_frst_res(num_frags)-1)=3
                  i_in_dom_fasta_tot(nres1+my_frst_res(num_frags)-1)=
     &              iam_in_pdb_chn(nnn,mmm,nres2)
                  my_nnn(nres1)=nnn
                  my_mmm(nres1)=mmm
                  my_nres2(nres1)=nres2
                  my_nres2_tru(nres1)=nres2_tru
                  nres_iden(num_chns_tot)=nres_iden(num_chns_tot)+1 ! v1.4
                endif
              elseif(streng0(21+i:21+i).eq.':') then
                if(iam_in_pdb_chn(nnn,mmm,nres2).ge.0) then
                  i_in_pdb_fasta(nres1)=2 ! similar match
                  i_in_pdb_fasta_tot(nres1+my_frst_res(num_frags)-1)=2
                  i_in_dom_fasta_tot(nres1+my_frst_res(num_frags)-1)=
     &              iam_in_pdb_chn(nnn,mmm,nres2)
                  my_nnn(nres1)=nnn
                  my_mmm(nres1)=mmm
                  my_nres2(nres1)=nres2
                  my_nres2_tru(nres1)=nres2_tru
                  nres_posi(num_chns_tot)=nres_posi(num_chns_tot)+1 ! v1.4
                endif
              elseif(streng0(21+i:21+i).eq.'.') then
                if(iam_in_pdb_chn(nnn,mmm,nres2).ge.0) then
                  i_in_pdb_fasta(nres1)=1 ! any kind of match
                  i_in_pdb_fasta_tot(nres1+my_frst_res(num_frags)-1)=1
                  i_in_dom_fasta_tot(nres1+my_frst_res(num_frags)-1)=
     &              iam_in_pdb_chn(nnn,mmm,nres2)
                  my_nnn(nres1)=nnn
                  my_mmm(nres1)=mmm
                  my_nres2(nres1)=nres2
                  my_nres2_tru(nres1)=nres2_tru
                  nres_alin(num_chns_tot)=nres_alin(num_chns_tot)+1 ! v1.4
                endif
              endif
              write(*,*)"ALIGN ",streng1(21+i:21+i),nres1,
     &                           streng2(21+i:21+i),nres2

C also store information for making an easy to read alignment file that
C will be written to falign2 subsequently

              ntot=ntot+1
              sequenceA(ntot)=streng1(21+i:21+i)
              sequenceM(ntot)=streng0(21+i:21+i)
              sequenceB(ntot)=streng2(21+i:21+i)

C v1.5 I tried removing the requirement that sequenceS,D,T only get written
C when they are aligned with the target but this was a bad idea - I
C think the reason I originally added the requirement was that it made
C sure that we didn't erronesouly write out * or @ entries for residues
C that are present in the target but absent in the template - if I
C remove the requirement then we continue to write out the previous *, @
C (and DSSP designations) for each residue until nres2 updates. The
C bottom line is that the original decision to add the "only if aligned"
C requirement was the correct one - hence it's reinstated below...

              if(nres2.ge.1) then
                if(iam_in_SS_elem(nnn,mmm,nres2).ge.1.and.
     &             streng2(21+i:21+i).ne.'-') then ! only if aligned
ccc             if(iam_in_SS_elem(nnn,mmm,nres2).ge.1) then
                  sequenceS(ntot)='*'

C v2.6 only set i_in_SS_elemnt & i_in_SS_elemnt_tot if nres1>0
 
                  if(nres1.ge.1) then
                    i_in_SS_elemnt(nres1)=1
                    i_in_SS_elemnt_tot(nres1+my_frst_res(num_frags)-1)=1
                  endif

                else
                  sequenceS(ntot)=' '
                endif
                if(my_DSSP_design(nnn,mmm,nres2).ne.' '.and. ! v1.4
     &             streng2(21+i:21+i).ne.'-') then ! only if aligned
ccc             if(my_DSSP_design(nnn,mmm,nres2).ne.' ') then
                  sequenceD(ntot)=my_DSSP_design(nnn,mmm,nres2)
                else
                  sequenceD(ntot)=' '
                endif
                if(iam_in_pdb_chn(nnn,mmm,nres2).ge.0.and.
     &             streng2(21+i:21+i).ne.'-') then ! only if aligned
ccc             if(iam_in_pdb_chn(nnn,mmm,nres2).eq.1) then
                  sequenceT(ntot)='@'
                else
                  sequenceT(ntot)=' '
                endif
              endif

C look for the beginning of deletions
C note that we discount cases where we haven't yet found a true res

              if(streng1(21+i:21+i).eq.'-'.and.
     &           streng2(21+i:21+i).ne.'-') then
                if(nres1.ge.1.and.indel.eq.0) then
                  idel_beg=nres1
                  indel=1
                endif
              endif

C look for end of deletions - note that we only mark residues as being
C in deletions when we have found both the beginning and the end...

              if(streng1(21+i:21+i).ne.'-'.and.
     &           indel.eq.1) then
                idel_end=nres1 ! this used to say nres1-1
                indel=0
                i_am_in_a_del(idel_beg)=1
                i_am_in_a_del(idel_end)=1
              endif

C look at the next character on the line

            enddo

C go back and read another block of lines

            goto 109

C come here when done with the alignment file

309         close(11)

C if use_SS_constraints and ifrst_time_thru then we make a note of
C len_pdb_chn_fin(nnn,mmm) - this can be used to truncate the template
C we do this by finding the max value of my_nres2 assigned to any of the
C nres1 residues...
C v1.4

            if(use_SS_constraints.and.ifrst_time_thru.eq.1) then
              my_nres2_max=0
              do nt=1,nres1
                my_nres2_max=max(my_nres2_max,my_nres2(nt))
              enddo
              if(my_nres2_max.gt.len_of_pdb_chn_fin(nnn,mmm)) then
                len_of_pdb_chn_fin(nnn,mmm)=my_nres2_max
              endif

C v1.6 - I think we need to override this truncation when it cuts into a
C secondary structure element... if so, then we keep incrementing until
C there is no SS element aligned in the template

C note that it may eventually be better to extend this so it's not only
C cases where the truncation cuts into a SS element - in some cases the
C unconstrained alignment might be so bad that it finishes right before
C an important secondary structure element - so perhaps we might want to
C find the next SS element and keep going till that one ends - the
C danger with that approach is that SS element might screw the alignment
C up completely because it might be one that we *don't* want!!
           
              if(iam_in_SS_elem(nnn,mmm,my_nres2_max).ge.1) then
308             my_nres2_max=my_nres2_max+1
ccc             write(*,*)'truncation override: for pdb/chn nnn,mmm'
                if(my_nres2_max.le.nres2) then
                  if(iam_in_SS_elem(nnn,mmm,my_nres2_max).eq.0) then
                    len_of_pdb_chn_fin(nnn,mmm)=my_nres2_max
                    goto 310
                  else
                    goto 308
                  endif
                endif
              endif
310           continue 

            endif 

C also in v1.4 we'll find how many res in the target occurred *before*
C and *after* the template sequence's alignment - at this point we
C already know how many were before (see nres_befr() above) - to find
C out how many were after we need to count back from the last target
C res, continually asking if it corresponded to a res in the template

            if(my_nres2_tru(nres1).ne.0) then
              nres_aftr(num_chns_tot)=-1
            else
              lt=1
              do nt=nres1-1,1,-1
                if(my_nres2_tru(nt).ne.0) then
                  nres_aftr(num_chns_tot)=lt
                  exit
                endif
                lt=lt+1
              enddo
            endif

C also in v1.4 we'll find how many res in the template occurred *before*
C and *after* the target sequence's alignment...

            nres_tmplte_befr(num_chns_tot)=my_nres2_tru(1)-1
            my_nres2_tru_max=0
            do n3=1,nres1
              my_nres2_tru_max=max(my_nres2_tru_max,my_nres2_tru(n3))
            enddo
            nres_tmplte_aftr(num_chns_tot)=nres2_tru-my_nres2_tru_max
            if(nres_tmplte_aftr(num_chns_tot).eq.0) 
     &         nres_tmplte_aftr(num_chns_tot)=-1 ! set to -1 if N/A
            write(*,*)'HELP ME ',mmm,nres1,nres2,
     &                        nres2_tru,my_nres2_tru_max

C first, we should check that the latest alignment had the correct
C number of residues in it - should match initial fasta number

            if(nres1.ne.num_res_in_seq(iii)) then
              write(*,*)'HANG ON ',nres1,num_res_in_seq(iii),
     &                   falign(iii,nnn,mmm)
              write(*,*)'first number should match second number above'
              write(*,*)
              stop
            endif

C second, we should be able to write out the alignment in a more easy to
C read format that also shows the SS elements and the resolved residues
C eventually I think the following information should also be included
C in the final AHEMODEL pdb file...
 
            flen=len(trim(falign(iii,nnn,mmm)))
            falign2(iii,nnn,mmm)(1:flen)=falign(iii,nnn,mmm)
            falign2(iii,nnn,mmm)(flen+1:flen+10)='.rewritten'
            flen=flen+10
            write(*,*)'check falign2 ',falign2(iii,nnn,mmm)(1:flen)

            open(unit=21,file=falign2(iii,nnn,mmm)(1:flen),
     &                   status='unknown')
            if(iii.eq.1) then
              write(25,*)
              write(25,*)
            endif

C v4.1 clarified following write so that we're clear that we mean input
C chain here (i.e. in order listed in input file)

            write(25,996)mmm,falign2(iii,nnn,mmm)(1:6)
996         format('final seq alignment for input-chain # ',i8,1x,a6)
C996         format('final sequence alignment for chain # ',i8)

C v3.3 at this stage we used to read PSIPRED prediction for this 
C sequence from Tim's flat files - but we now remove all this and go
C straight for the .ss2 file

C v3.0 here is the original path to Tim's results

c           sys_string1( 1:33)='/home/tcollingsworth/PSI-RESULTS_'
c           sys_string1(34:44)='FOR_VENKAT/'
c           sys_string1(45:50)=falign2(iii,nnn,mmm)(1:6)
c           sys_string1(51:51)='/'
c           sys_string1(52:57)=falign2(iii,nnn,mmm)(1:6)
c           sys_string1(58:61)='.ss2'

C here is the new path - note that for backwards compatibility we
C require that each .ss2 file be in a directory with the same filename
C stem as the .ss2 file...

C v3.9 make sure that ss2 file name is correct when the sequence name is
C not six characters - nname_seq(n) is the length of the filename stem 

C v4.7 it is now assumed that ss2_path is the full path and file name of
C the .ss2 file - so there should be no need to add any crap to the end
C of it - hence the commented out lines that add a / and filename below

C 2022 v1.0 THAT MEANS IT ONLY WORKS FOR ONE CHAIN FFS - MUCH BETTER TO
C JUST DROP THE PATH ENTIRELY AND REQUIRE USERS TO HAVE P12345.SS2 TO
C ACCOMPANY P12345.FASTA IF THEY WANT TO USE SS2 FILES...

            sys_string1(1:6)=falign2(iii,nnn,mmm)(1:6)
            sys_string1(7:10)='.ss2'

C v2.6 - if we set the uniprot fasta to PXXXXX then we will not go
C looking for a .ss2 file...

C v3.7 - UPDATE we now use the inquire intrinsic to avoid crashes when
C the ss2 file doesn't exist...

C v3.9 also use this to set have_ss2_file to true - have_ss2_file=true
C if there is so much as one ss2 file in existence...

            ss2_file_exists=.false.
            inquire(file=sys_string1(1:10),exist=ss2_file_exists)
            if(ss2_file_exists) have_ss2_file=.true.

            if(ss2_file_exists) then
              write(*,*)
              write(*,*)'Now reading .ss2 file ',
     &                   sys_string1(1:10)
              write(*,*)
              open(unit=29,file=sys_string1(1:10),status='unknown')

C first, we will try to find whether the .ss2 sequence is a match for
C the current sequence - so we read the first 10 AAs from the .ss2 -
C note that we allow for missing or extra residue in the .ss2 value
C through use of "noffset" - I think this was originally implemented for
C signal peptides (present in the .ss2 file but not in the .fasta???)

              read(29,'(1x)')
              read(29,'(1x)')
              lll=0
4775          read(29,'(a30)',end=4776)string30
              lll=lll+1
              read(string30,'(5x,a)')char1
              seqss2(lll:lll)=char1
              if(lll.eq.10) goto 4776
              goto 4775
4776          rewind(29)

              noffset=-999
              do i22=1,num_res_in_seq(iii)-10
                imatch=1
                do j22=1,10
                  if(seq_of_fasta(iii)(i22+j22-1:i22+j22-1).ne.
     &               seqss2(j22:j22)) imatch=0
                enddo
                if(imatch.eq.1) then
                  noffset=i22-1
                  exit
                endif
              enddo
              if(noffset.eq.-999) then
                write(*,*)'sequence mismatch in .ss2 file '
                stop
              else
                write(*,*)'ss2 file offset = ',noffset
              endif

              read(29,'(1x)')
              read(29,'(1x)')
              lll=0
5775          read(29,'(a30)',end=5776)string30
              lll=lll+1
              ll2=lll+noffset

C v4.2 when running nascent protein chains the 4.1 code crashes when the
C ss2 file is full-length - we want the code to continue even if the
C .ss2 file is *longer* than the input .fasta so let's see if the
C following line catches those problems...

              if(ll2.gt.ntot) then
                write(*,*)
                write(*,*)'WARNING - .ss2 file is longer than .fasta'
                write(*,*)'IGNORE this if the chain is a nascent one'
                write(*,*)
                goto 5776
              endif

              read(string30,'(5x,a,1x,a,5x,a,6x,a,6x,a)')char1,
     &                                          ajunk9,b1,b2,b3

C v3.5 change the (incorrect) code that checks whether the .ss2 sequence
C matches that for the .fasta - we used to compare char1 with
C sequenceA(ll2) but the latter comes from an alignment file and can
C contain "-" and " " entries - so now we do something really clunky and
C we just read through sequenceA till we get to the corresponding
C sequence position and then check - this is appallingly clunky code but
C it should now work: ll2 is the actual residue number, ll3 is the
C position in the alignment...

              ll4=0
              do ll3=1,ntot
                if(sequenceA(ll3).ne.' '.and.
     &             sequenceA(ll3).ne.'-') ll4=ll4+1
                if(ll4.eq.ll2) exit ! we now know ll3, so use it below
              enddo

C note that we only need "ll3" for the sequence comparison test - for
C all subsequent statements it's correct to use "ll2" as this is the
C actual residue number, not the position in the sequence alignment...
              
              if(sequenceA(ll3).ne.char1) then
                write(*,*)
                write(*,*)'sequence mismatch in .ss2 file '
                write(*,*)'iii = ',iii
                write(*,*)'nnn = ',nnn
                write(*,*)'lll = ',lll
                write(*,*)'seq = ',sequenceA(ll3)
                write(*,*)'cha = ',char1
                write(*,*)'quitting :('
                write(*,*)
                stop
              endif
              my_PSIPRED_prd(iii,nnn,ll2)=ajunk9 ! correct to use ll2
              if(ajunk9.eq.'C') then             ! here and below...
                my_PSIPRED_cnf(iii,nnn,ll2)=b1
              elseif(ajunk9.eq.'H') then
                my_PSIPRED_cnf(iii,nnn,ll2)=b2
              elseif(ajunk9.eq.'E') then
                my_PSIPRED_cnf(iii,nnn,ll2)=b3
              endif
              read(string30,'(9x,3f7.3)')f1,f2,f3
              my_PSIPRED_vl1(iii,nnn,ll2)=f1
              my_PSIPRED_vl2(iii,nnn,ll2)=f2
              my_PSIPRED_vl3(iii,nnn,ll2)=f3
              goto 5775
5776          close(29)
            else
              write(*,*)
              write(*,*)'WARNING: ',sys_string1(1:10),' does not exist'
              write(*,*)
            endif

            nres_tmp=0
            do lll=1,ntot
              if(sequenceA(lll).ne.' '.and.
     &           sequenceA(lll).ne.'-') then
                nres_tmp=nres_tmp+1
                sequenceP(lll:lll)=my_PSIPRED_prd(iii,nnn,nres_tmp)
                sequenceC(lll:lll)=my_PSIPRED_cnf(iii,nnn,nres_tmp)

C if this is the first match for this sequence then we write out the
C DSSP and PSIPRED values for future analysis...

                if(ididthisseq.eq.1) then
                  write(93,9101)iii,nres_tmp,sequenceD(lll:lll),
     &                          my_PSIPRED_prd(iii,nnn,nres_tmp),
     &                          my_PSIPRED_vl1(iii,nnn,nres_tmp),
     &                          my_PSIPRED_vl2(iii,nnn,nres_tmp),
     &                          my_PSIPRED_vl3(iii,nnn,nres_tmp)
9101              format('seq# ',i4,' res# ',i6,
     &                   '   DSSP = ',a,'   PSIPRED = ',a,
     &         '     PSIPRED vals: C = ',f7.3,' H = ',f7.3,' E = ',f7.3)
                endif

C v2.9 - for PSIPRED-predicted helices we take the helical probability,
C multiply by 2 and subtract 1 - this ensures that low confidence
C predictions (prob = 0.5) end up with occupancies near zero - do
C something similar with PSIPRED-predicted sheets but make negative...

                nres_final=nres_final+1
                if(my_PSIPRED_vl2(iii,nnn,nres_tmp).ge.
     &             my_PSIPRED_vl1(iii,nnn,nres_tmp).and.
     &             my_PSIPRED_vl2(iii,nnn,nres_tmp).ge.
     &             my_PSIPRED_vl3(iii,nnn,nres_tmp)) then
                  occ_psipred(nres_final)=
     &              my_PSIPRED_vl2(iii,nnn,nres_tmp)*2.0-1.0
                endif
                if(my_PSIPRED_vl3(iii,nnn,nres_tmp).ge.
     &             my_PSIPRED_vl1(iii,nnn,nres_tmp).and.
     &             my_PSIPRED_vl3(iii,nnn,nres_tmp).ge.
     &             my_PSIPRED_vl2(iii,nnn,nres_tmp)) then
                  occ_psipred(nres_final)=
     &             -my_PSIPRED_vl3(iii,nnn,nres_tmp)*2.0+1.0
                endif
     
              endif
            enddo
 
            write(21,*)
            write(25,*)
            nrounds=int(ntot/100)+1
            do n=1,nrounds
              ibeg=(n-1)*100+1
              iend=n*100
              if(n.eq.nrounds) then
                do m=ntot+1,iend
                  sequenceA(m)=' ' ! sequence of target
                  sequenceM(m)=' ' ! shows identity/similarity
                  sequenceB(m)=' ' ! sequence of template
                  sequenceS(m)=' ' ! shows whether in SS element
                  sequenceD(m)=' ' ! shows DSSP designation
                  sequenceP(m)=' ' ! shows PSIPRED prediction
                  sequenceC(m)=' ' ! shows PSIPRED confidence
                  sequenceT(m)=' ' ! shows whether in pdb file
                  sequenceU(m)=' ' ! shows total res#
                enddo
              endif
              write(21,6611)(sequenceA(ibeg:iend))
              write(21,6612)(sequenceM(ibeg:iend))
              write(21,6613)(sequenceB(ibeg:iend))
              write(21,6614)(sequenceS(ibeg:iend))
              write(21,6615)(sequenceD(ibeg:iend))
              write(21,6616)(sequenceP(ibeg:iend))
              write(21,6617)(sequenceC(ibeg:iend))
              write(21,6618)(sequenceT(ibeg:iend))
c             write(21,6619)(sequenceU(ibeg:iend))
              write(21,*)
              write(21,*)
              write(25,6611)(sequenceA(ibeg:iend))
              write(25,6612)(sequenceM(ibeg:iend))
              write(25,6613)(sequenceB(ibeg:iend))
              write(25,6614)(sequenceS(ibeg:iend))
              write(25,6615)(sequenceD(ibeg:iend))
              write(25,6616)(sequenceP(ibeg:iend))
              write(25,6617)(sequenceC(ibeg:iend))
              write(25,6618)(sequenceT(ibeg:iend))
c             write(25,6619)(sequenceU(ibeg:iend))
              write(25,*)
              write(25,*)
6601          format(100a)
6611          format(100a,'  <-- target sequence')
6612          format(100a,'  <-- alignment match')
6613          format(100a,'  <-- template sequence')
6614          format(100a,'  <-- SS used to align')
6615          format(100a,'  <-- template DSSP SS ')
6616          format(100a,'  <-- target PSIPRED SS')
6617          format(100a,'  <-- target PSIPRED conf')
6618          format(100a,'  <-- res in template')
6619          format(100a,'  <-- res in template')
6602          format(100i1)
            enddo
            close(21)

C new for v1.3 - also write out in rosetta format...

            falign2(iii,nnn,mmm)(flen+1:flen+8)='.rosetta'
            flen=flen+8
            write(*,*)'check falign2 ',falign2(iii,nnn,mmm)(1:flen)

            open(unit=21,file=falign2(iii,nnn,mmm)(1:flen),
     &                   status='unknown')
            write(21,'(">",a6)')falign2(iii,nnn,mmm)(1:6)
            nrounds=int(ntot/100)+1
            do n=1,nrounds
              ibeg=(n-1)*100+1
              iend=n*100
              if(n.eq.nrounds) then
                do m=ntot+1,iend
                  sequenceA(m)=' ' ! sequence of target
                enddo
              endif
              write(21,6601)(sequenceA(ibeg:iend))
            enddo
            write(21,'(">TEMPLATE.pdb")')
            do n=1,nrounds
              ibeg=(n-1)*100+1
              iend=n*100
              if(n.eq.nrounds) then
                do m=ntot+1,iend
                  sequenceB(m)=' ' ! sequence of template
                enddo
              endif

C we need to make sure to write in hyphens for missing residues

              do m=ibeg,iend
                if(sequenceB(m).ne.' ') then
                  if(sequenceT(m).ne.'@') sequenceB(m)='-'
                endif
              enddo

              write(21,6601)(sequenceB(ibeg:iend))
            enddo
            close(21)

C at this point we have attempted to match all residues in our sequence
C with the pdb/chain combo - we should now be able to write out a pdb
C file containing the matched residues...

C now let's figure out where the tails and loops are...

            len_nter=0
            do i=1,nres1 
              if(i_in_pdb_fasta(i).eq.0) then
                len_nter=len_nter+1
              else
                exit
              endif
            enddo
            write(*,*)'length of Nter = ',len_nter

            len_cter=0
            do i=nres1,1,-1
              if(i_in_pdb_fasta(i).eq.0) then
                len_cter=len_cter+1
              else
                exit
              endif
            enddo
            write(*,*)'length of Cter = ',len_cter

C here we should check the template pdb to make sure there aren't going
C to be big gaps between aligned residues that are adjacent in target
C start by reading the N and C atoms of each residue in template pdb
C UPDATE - this code currently commented out - probably unnecessary...

            nres2=0
            open(unit=11,file=fname_pdb_chn(nnn,mmm),status='unknown')
151         read(11,21,end=152)strong
            read(strong(23:26),'(i4)')nres2
            if(strong(13:16).eq.' N  ') 
     &         read(strong,153)xntm(nres2),yntm(nres2),zntm(nres2)
            if(strong(13:16).eq.' C  ') 
     &         read(strong,153)xctm(nres2),yctm(nres2),zctm(nres2)
            goto 151
153         format(30x,3f8.3)
152         close(11)
 
            num_loops=0
c           do i=1,nres1-1
c             if(i_in_pdb_fasta(i).ge.1.and.
c    &           i_in_pdb_fasta(i+1).ge.1) then
c               itmp1=my_nres2(i)
c               itmp2=my_nres2(i+1)
c               dist=sqrt(
c    &                   (xctm(itmp1)-xntm(itmp2))**2+
c    &                   (yctm(itmp1)-yntm(itmp2))**2+
c    &                   (zctm(itmp1)-zntm(itmp2))**2)
c               if(dist.gt.2.0) then ! it's a peptide bond so 2A
c                 num_loops=num_loops+1
c                 ibeg_loop(num_loops)=i
c                 iend_loop(num_loops)=i+1
c                 write(*,*)
c                 write(*,*)'WARNING!!! disconnect found in template'
c                 write(*,*)'in file ',fname_pdb_chn(nnn,mmm)
c                 write(*,*)'between tmplte residues ',itmp1,' & ',itmp2
c                 write(*,*)'between target residues ',i,' & ',i+1
c                 write(*,*)'distance = ',dist
c       write(*,*)'itmp1 coord ',xctm(itmp1),yctm(itmp1),zctm(itmp1)
c       write(*,*)'itmp2 coord ',xntm(itmp2),yntm(itmp2),zntm(itmp2)
c                 write(*,*)'we will add a loop here...',num_loops
c                 write(*,*)
c               endif
c             endif
c           enddo

C here's where we find the insertions: we just look for residues that
C match to a position that has no coords (use i_in_pdb_fasta for this)

            in_loop=0
            do i=len_nter+1,nres1-len_cter
              if(i_in_pdb_fasta(i).eq.0) then
                if(in_loop.eq.0) then
                  in_loop=1
                  num_loops=num_loops+1
                  write(*,*)'new loop insertion@ ',num_loops,i
                  ibeg_loop(num_loops)=i
                  iend_loop(num_loops)=i
                else
                  iend_loop(num_loops)=i
                endif
              else
                in_loop=0
              endif
            enddo

            write(*,*)'number of insertion loops = ',num_loops

C we can use basically the same code to find deletions and treat them as
C regular loops for now - not sure if this will work but worth a shot
C UPDATE - I think this code works okay - if we use this code then I
C think we basically don't need to check connection lengths (see above)

            in_loop=0
            do i=len_nter+1,nres1-len_cter
              if(i_am_in_a_del(i).eq.1) then
                write(*,*)'in a del ',i,in_loop
                if(in_loop.eq.0) then
                  in_loop=1
                  num_loops=num_loops+1
                write(*,*)'new loop deletion ',num_loops
                  ibeg_loop(num_loops)=i-1
                  iend_loop(num_loops)=i
                  write(*,*)'found deletion starting at ',i-1
                  write(*,*)'deletion assumed to end at ',i
                else
                  iend_loop(num_loops)=i
                  write(*,*)'deletion now ends at ',i
                endif
              else
                in_loop=0
              endif
            enddo

            write(*,*)'FOR ALIGNMENT ',falign(iii,nnn,mmm)
            write(*,*)'BEFORE TWEAKING LOOPS'
            do i=1,num_loops

C v2.6 oh great... - let's fix messed up loops again...
         
              if(ibeg_loop(i).lt.1)     ibeg_loop(i)=1
              if(iend_loop(i).gt.nres1) iend_loop(i)=nres1

              write(*,*)'LOOP ',i,ibeg_loop(i),iend_loop(i),
     &          seq_of_fasta(iii)(ibeg_loop(i):iend_loop(i))
            enddo

C now we may need to amalgamate loops... - we do this in iterative
C fashion - just doing pairs at a time until no further changes

607         num_loops_final=0
            idone=0
            do i=1,num_loops
              if(idone(i).eq.1) cycle ! skips 2nd of the merged loops
              iaccept=1
              do j=i+1,num_loops

C we just need to look for overlapping loops - it turns out that it's
C easier to identify those that cannot possibly overlap and then treat
C all others as overlapping... note that once we have one overlap for
C the current loop (i), we don't try to do any more merging in this
C round - we don't try to fix everything in one go (too dangerous)

                ioverlap=1
                if((iend_loop(j).lt.ibeg_loop(i)-2)) ioverlap=0
                if((ibeg_loop(j).gt.iend_loop(i)+2)) ioverlap=0
                if(ioverlap.eq.0) cycle ! skip if no possible overlap
                iaccept=0
                idone(j)=1 ! don't do anything further with j
                num_loops_final=num_loops_final+1
                write(*,*)'merging two loops ',i,' & ',j
                if(ibeg_loop(j).le.ibeg_loop(i)) then
                  ibeg_loop_final(num_loops_final)=ibeg_loop(j)
                else
                  ibeg_loop_final(num_loops_final)=ibeg_loop(i)
                endif
                if(iend_loop(j).ge.iend_loop(i)) then
                  iend_loop_final(num_loops_final)=iend_loop(j)
                else
                  iend_loop_final(num_loops_final)=iend_loop(i)
                endif
                if(iaccept.eq.0) exit ! exit as soon as we have a merge
                                      ! this should fix cases where a
                                      ! loop merges with 2 or more...
              enddo
              if(iaccept.eq.1) then
                num_loops_final=num_loops_final+1
                write(*,*)'accepting loop ',i,' unchanged'
                ibeg_loop_final(num_loops_final)=ibeg_loop(i)
                iend_loop_final(num_loops_final)=iend_loop(i)
              endif
            enddo
            if(num_loops_final.ne.num_loops) then
              write(*,*)'rechecking merges ',num_loops_final,num_loops
              do i=1,num_loops_final
                ibeg_loop(i)=ibeg_loop_final(i)
                iend_loop(i)=iend_loop_final(i)
              enddo
              num_loops=num_loops_final
              goto 607
            endif

            write(*,*)'AFTER MERGING LOOPS THE FIRST TIME'
            do i=1,num_loops
              write(*,*)'LOOP ',i,ibeg_loop(i),iend_loop(i),
     &          seq_of_fasta(iii)(ibeg_loop(i):iend_loop(i))
            enddo

C before going on, immediately tweak insertions in cases where the
C residues before and after the loop are structured and contiguous in
C the final sequence... for now we'll just add both residues to either
C side to the loop - I guess we could do something fancier (see
C commented out part) later but this should be okay for now...

            do i=1,num_loops
              ibeg0=ibeg_loop(i)-1
              iend0=iend_loop(i)+1

C v2.6 need to adjust these numbers if they take use too far - note that
C they should eventually get changed to tails rather than loops...

              if(ibeg0.lt.1) then
                ibeg0=1
                iend0=iend0+1
              endif
              if(iend0.gt.nres1) then
                ibeg0=ibeg0-1
                iend0=nres1
              endif

              if(i_in_pdb_fasta(ibeg0).ne.0.and.
     &           i_in_pdb_fasta(iend0).ne.0.and.
     &           my_nres2(iend0).eq.my_nres2(ibeg0)+1) then
                write(*,*)
                write(*,*)'loop # ',i,' is expecting too much! '
                write(*,*)'cannot fit in a loop between adjacent res'
                write(*,*)'so increase the loop on either side'
                write(*,*)
                ibeg_loop(i)=ibeg0
                iend_loop(i)=iend0
c               if(i_in_pdb_fasta(ibeg_loop(i)).eq.
c    &             i_in_pdb_fasta(iend_loop(i))) then
c                 call random_number(u1)
c                 if(u1.le.0.5) then
c                   ibeg_loop(i)=ibeg_loop(i)-1
c                 else
c                   iend_loop(i)=iend_loop(i)+1
c                 endif
c               elseif(i_in_pdb_fasta(ibeg_loop(i)).gt.
c    &                 i_in_pdb_fasta(iend_loop(i))) then
c                 iend_loop(i)=iend_loop(i)+1
c               elseif(i_in_pdb_fasta(ibeg_loop(i)).lt.
c    &                 i_in_pdb_fasta(iend_loop(i))) then
c                 ibeg_loop(i)=ibeg_loop(i)-1
c               endif
              endif
            enddo

C first we need to make sure that loops are at least min_loop_size residues
C note that we compare the beg and end residues of the loop and extend
C beyond the one that has the lowest level of match with the template -
C if they have equal priority we pick a random number to choose...

C v1.9 before doing this, we should record the current states of the 
C loops - we can then compare these with previous chains in order to
C make sure that all loops that are common to other chains are treated
C in exactly the same way...

            do i=1,num_loops
              ibeg_ori_loop(i)=ibeg_loop(i)
              iend_ori_loop(i)=iend_loop(i)
            enddo

            do i=1,num_loops
60606         continue

C v2.6 immediately skip out of here if the loop is at the end of
C sequence - we don't want to do silly stuff like expand the loop
C for something that's destined to be a tail anyway...

              if(ibeg_loop(i).eq.1) cycle
              if(iend_loop(i).eq.nres1) cycle

C v4.1 here we try to add a fix for cases where we have a single placed
C residue at either the N-terminus or the C-terminus - I'm thinking of
C cases like making lipidated proteins where we'll have a single CYS
C residues at the N-terminus - we'll copy the lines from above but apply
C them to the 2nd and penultimate residues respectively...

              if(rebuild_auto) then
                if(ibeg_loop(i).eq.2) cycle
                if(iend_loop(i).eq.nres1-1) cycle
              endif

              idiff=iend_loop(i)-ibeg_loop(i)+1
              if(idiff.lt.min_loop_size) then

C v1.9
C if we failed to build the loop we need to increase its size in the
C hope of building it successfully - here we make use of
C i_in_pdb_fasta_tot which helps us know which of the positions to
C extend the loop in...

                if(i_in_SS_elemnt(ibeg_loop(i)).gt.
     &             i_in_SS_elemnt(iend_loop(i))) then
                  iend_loop(i)=iend_loop(i)+1
                  write(*,*)'N-ter loop res is in SS elemnt so survives'
                  goto 60606
                elseif(i_in_SS_elemnt(ibeg_loop(i)).lt.
     &                 i_in_SS_elemnt(iend_loop(i))) then
                  ibeg_loop(i)=ibeg_loop(i)-1
                  write(*,*)'C-ter loop res is in SS elemnt so survives'
                  goto 60606
                endif

C if both the beg and end were already in SS elements then we'll need to
C use sequence conservation to break the tie...

                if(i_in_pdb_fasta(ibeg_loop(i)).eq.
     &             i_in_pdb_fasta(iend_loop(i))) then
                  call random_number(u1)
                  if(u1.le.0.5) then
                    ibeg_loop(i)=ibeg_loop(i)-1 
                  else
                    iend_loop(i)=iend_loop(i)+1 
                  endif
                elseif(i_in_pdb_fasta(ibeg_loop(i)).gt.
     &                 i_in_pdb_fasta(iend_loop(i))) then
                  iend_loop(i)=iend_loop(i)+1 
                elseif(i_in_pdb_fasta(ibeg_loop(i)).lt.
     &                 i_in_pdb_fasta(iend_loop(i))) then
                  ibeg_loop(i)=ibeg_loop(i)-1 
                endif
                goto 60606
              endif

C once we have made sure that the loop is of at least the min. size, we
C make sure that the loop didn't go beyond the limits of the chain if it
C got lengthened to make the minimum # of residues... note that we
C remove such loops later and make them Nter and Cter tails...

C v2.6 no need for these lines now since we now "cycle" above...

ccc           if(ibeg_loop(i).lt.1) ibeg_loop(i)=1
ccc           if(iend_loop(i).gt.nres1) iend_loop(i)=nres1

            enddo

            write(*,*)'AFTER  TWEAKING LOOPS'
            do i=1,num_loops
              write(*,*)'LOOP ',i,ibeg_loop(i),iend_loop(i),
     &          seq_of_fasta(iii)(ibeg_loop(i):iend_loop(i))
            enddo

C now we may again need to amalgamate loops... - we do this in iterative
C fashion - just doing pairs at a time until no further changes

707         num_loops_final=0
            idone=0
            do i=1,num_loops
              if(idone(i).eq.1) cycle ! skips 2nd of the merged loops
              iaccept=1
              do j=i+1,num_loops

C we just need to look for overlapping loops - it turns out that it's
C easier to identify those that cannot possibly overlap and then treat
C all others as overlapping... note that once we have one overlap for
C the current loop (i), we don't try to do any more merging in this
C round - we don't try to fix everything in one go (too dangerous)

                ioverlap=1
                if((iend_loop(j).lt.ibeg_loop(i)-2)) ioverlap=0
                if((ibeg_loop(j).gt.iend_loop(i)+2)) ioverlap=0
                if(ioverlap.eq.0) cycle ! skip if no possible overlap
                iaccept=0
                idone(j)=1 ! don't do anything further with j
                num_loops_final=num_loops_final+1
                write(*,*)'merging two loops ',i,' & ',j
                if(ibeg_loop(j).le.ibeg_loop(i)) then
                  ibeg_loop_final(num_loops_final)=ibeg_loop(j)
                else
                  ibeg_loop_final(num_loops_final)=ibeg_loop(i)
                endif
                if(iend_loop(j).ge.iend_loop(i)) then
                  iend_loop_final(num_loops_final)=iend_loop(j)
                else
                  iend_loop_final(num_loops_final)=iend_loop(i)
                endif
                if(iaccept.eq.0) exit ! exit as soon as we have a merge
                                      ! this should fix cases where a
                                      ! loop merges with 2 or more...
              enddo
              if(iaccept.eq.1) then
                num_loops_final=num_loops_final+1
                write(*,*)'accepting loop ',i,' unchanged'
                ibeg_loop_final(num_loops_final)=ibeg_loop(i)
                iend_loop_final(num_loops_final)=iend_loop(i)
              endif
            enddo
            if(num_loops_final.ne.num_loops) then
              write(*,*)'rechecking merges ',num_loops_final,num_loops
              do i=1,num_loops_final
                ibeg_loop(i)=ibeg_loop_final(i)
                iend_loop(i)=iend_loop_final(i)
              enddo
              num_loops=num_loops_final
              goto 707
            endif

            write(*,*)'AFTER MERGING LOOPS FOR THE SECOND TIME'
            do i=1,num_loops
              write(*,*)'LOOP ',i,ibeg_loop(i),iend_loop(i),
     &          seq_of_fasta(iii)(ibeg_loop(i):iend_loop(i))
            enddo

C now we can express all loops & tails in terms of total residues...
C note that we also store chain ID for loopy etc (ichn_tot_nters...)
C note that we also store chain-local numbers so that we can easily tell
C when a newly-extended loop encroaches on the Nter or Cter of the same
C chain later...

C first, store info for all residues including non-loop, non-tail res's

            num_tot_res_old=num_tot_res
            num_loc_res=num_res_in_seq(iii)
            do i=1,num_res_in_seq(iii)
              num_tot_res=num_tot_res+1

C v2.3 we need a 2D array to convert between chainID and res# that might
C be found in a pdb file and the actual total residue number...

              do j=1,chainIDmax
                if(cur_chainID(num_frags).eq.blah(j:j)) then
                  chainlocal_2_totres(j,i+cur_offset(num_frags))=
     &              num_tot_res
                endif
              enddo

              seq_tot_res(num_tot_res:num_tot_res)=
     &          seq_of_fasta(iii)(i:i)
c v2.1        ijk=mod(num_frags,chainIDmax)
c v2.1        if(ijk.eq.0) ijk=chainIDmax
c v2.1        chn_tot_res(num_tot_res:num_tot_res)=blah(ijk:ijk)
              chn_tot_res(num_tot_res:num_tot_res)=
     &          fin_chainID(num_frags)
              res_tot_res(num_tot_res)=i
            enddo

C next, store info for all residues in the N-terminal tail if present

            if(len_nter.gt.0) then
              num_tot_nters=num_tot_nters+1
              ifrg_tot_nters(num_tot_nters)=num_frags
              ibeg_tot_nters(num_tot_nters)=num_tot_res_old+1 
              iend_tot_nters(num_tot_nters)=num_tot_res_old+len_nter
              ibeg_loc_nters(num_tot_nters)=1 
              iend_loc_nters(num_tot_nters)=len_nter
c v2.1        ijk=mod(num_frags,chainIDmax)
c v2.1        if(ijk.eq.0) ijk=chainIDmax
c v2.1        ichn_tot_nters(num_tot_nters)=blah(ijk:ijk)
              ichn_tot_nters(num_tot_nters)=cur_chainID(num_frags)
            endif

C now, store info for all residues in the C-terminal tail if present

            if(len_cter.gt.0) then
              num_tot_cters=num_tot_cters+1
              ifrg_tot_cters(num_tot_cters)=num_frags
              ibeg_tot_cters(num_tot_cters)=num_tot_res-len_cter+1
              iend_tot_cters(num_tot_cters)=num_tot_res
              ibeg_loc_cters(num_tot_cters)=num_loc_res-len_cter+1
              iend_loc_cters(num_tot_cters)=num_loc_res
c v2.1        ijk=mod(num_frags,chainIDmax)
c v2.1        if(ijk.eq.0) ijk=chainIDmax
c v2.1        ichn_tot_cters(num_tot_cters)=blah(ijk:ijk)
              ichn_tot_cters(num_tot_cters)=cur_chainID(num_frags)

C in order that OXT atoms be placed for residues that are at the actual
C C-terminus of each chain we need to store them also - note that this
C is only an issue when the input pdb file *already contains* the
C C-terminal residue - if it doesn't the code will build it in properly

            elseif(len_cter.eq.0) then
              num_act_cters=num_act_cters+1
              ipos_act_cters(num_act_cters)=num_tot_res
            endif

C finally, store info for all residues in loops if present - note that
C here we need to check if we've ever dealt with a chain of the same
C sequence before - if we have, then we want to make the same decisions
C about loop expansions and mergers that we made with the first chain...

            do i=1,num_loops
              num_tot_loops=num_tot_loops+1
              ifrg_tot_loops(num_tot_loops)=num_frags
              ibeg_tot_loops(num_tot_loops)=num_tot_res_old+
     &                                      ibeg_loop(i)
              iend_tot_loops(num_tot_loops)=num_tot_res_old+
     &                                      iend_loop(i)
              ibeg_loc_loops(num_tot_loops)=ibeg_loop(i)
              iend_loc_loops(num_tot_loops)=iend_loop(i)
              ibeg_ori_loops(num_tot_loops)=ibeg_ori_loop(i)
              iend_ori_loops(num_tot_loops)=iend_ori_loop(i)
c v2.1        ijk=mod(num_frags,chainIDmax)
c v2.1        if(ijk.eq.0) ijk=chainIDmax
c v2.1        ichn_tot_loops(num_tot_loops)=blah(ijk:ijk)
              ichn_tot_loops(num_tot_loops)=cur_chainID(num_frags)

C loop over all prior loops

              do j=1,num_tot_loops-1

C if they're from different types of chains then cycle

                if(my_frags_seq_num(ifrg_tot_loops(j)).ne.
     &             my_frags_seq_num(ifrg_tot_loops(num_tot_loops))) then
                  cycle
                endif

C otherwise check to see if their original loop assignments were
C identical with those of the current one...

                if(ibeg_ori_loops(j).eq.
     &             ibeg_ori_loops(num_tot_loops).and.
     &             iend_ori_loops(j).eq.
     &             iend_ori_loops(num_tot_loops)) then

C if so, then we need to correct ibeg_tot_loops & iend_tot_loops
C according to the difference of what was set for the earlier chain's
C version of this loop and what was set above for this loop

                  write(*,*)'overriding loop assignments for chain# ',
     &            num_frags,' loop# ',i,' to match identical chain# ',j
                  ibeg_tot_loops(num_tot_loops)=
     &            ibeg_tot_loops(num_tot_loops)+ibeg_loc_loops(j)-
     &                            ibeg_loc_loops(num_tot_loops)
                  iend_tot_loops(num_tot_loops)=
     &            iend_tot_loops(num_tot_loops)+iend_loc_loops(j)-
     &                            iend_loc_loops(num_tot_loops)

C we also need to set ibeg_loc_loops & iend_loc_loops so that they match
C those used for the earlier chain's version of this loop

                  ibeg_loc_loops(num_tot_loops)=ibeg_loc_loops(j)
                  iend_loc_loops(num_tot_loops)=iend_loc_loops(j)

C having found a match with an earlier chain's loop we should exit this
C do loop - shouldn't matter but no point taking a chance on this

                  exit
                endif
              enddo
            enddo

C at this point we should be able to start writing out to a grand pdb
C file and a grand fasta file in scwrl format so that we can run scwrl
C on the thing for all aligned residues...

C start by reading in the aligned pdb file...

C NOTE THAT WE ARE USING THE N ATOM HERE TO DETERMINE WHEN A NEW RESIDUE
C STARTS BUT WE'D LIKE TO USE CA - WE NEED TO SORT THIS OUT - PROBABLY
C THE BEST THING TO DO WILL BE TO DO A PROPER FILTERING OF ALL SHIT
C RESIDUES IN EACH OF THE PDB FILES BEFORE DOING ANYTHING ELSE...
C IT LOOKS LIKE SCWRL CAN HANDLE INCOMPLETE SIDECHAINS BUT NOT
C INCOMPLETE BACKBONE - SO WE DO NEED TO DO THIS AS FIRST THING

C we keep track of any rigid domains that are already assigned in the
C input PDB file by reading the occupancy field

            num_domains(num_frags)=0

            nres2=0
            open(unit=11,file=fname_pdb_chn(nnn,mmm),status='unknown')
47          read(11,21,end=48)strong

C read the occupancy field to keep track of num_domains

            read(strong,'(54x,f6.2)')occupancy
            iocc=int(occupancy)
            if(iocc.gt.num_domains(num_frags)) then
              num_domains(num_frags)=num_domains(num_frags)+1
            endif

C read the beta-factor field...

            read(strong,'(60x,f6.2)')beta
 
C override the chain ID so that scwrl doesn't reorder atoms later
C v2.1 the above comment is a bit cryptic but suggests that scwrl might
C reorder chains if read in as B then A etc - that's okay since we don't
C have that situation in our typical templates (can't rule it out
C though) - but it's also possible that if we recycle a chainID then we
C might have both chains with that chainID written out one after the
C other, which would be totally crap :(
C UPDATE - we've established that this is indeed the case: all chainA
C entries, for example, get lumped together - and in residue number
C order so god knows what would happen if the same residue number and
C chainID got used twice - scwrl would probably go crazy in writeout

C so, the final answer is to use the (code-local) cur_chainID from above

c v2.1      ijk=mod(num_frags,chainIDmax)
c v2.1      if(ijk.eq.0) ijk=chainIDmax
c v2.1      strong(22:22)=blah(ijk:ijk)
            strong(22:22)=cur_chainID(num_frags)

C note that nres_scwrl is the residue counter in the grand scwrl pdb file
C we need to make sure to pass loops and termini to this list also for
C subsequent use of loopy and our own tail generator...

            read(strong(23:26),'(i4)')nres2

C we need to compare all residues of the current sequence with the
C current residue found in this pdb file - if it's one that's a perfect
C match then we write out its backbone and sidechain and we set its
C fasta entry to be lower-case so that scwrl does nothing with it later
C - if it's an imperfect match then we just keep its backbone and ask
C scwrl to build it in later...

            do i=1,nres1

C do the following if the matching res in the pdb is *not similar*:

              if(i_in_pdb_fasta(i).eq.1.and.
     &           nres2.eq.my_nres2(i)) then
                if(strong(13:16).eq.' N  ') then
                  nres_scwrl=nres_scwrl+1
                  my_scwrl_res(num_tot_res_old+i)=nres_scwrl
                  my_tot_res(nres_scwrl)=num_tot_res_old+i

                  i_am_original(my_tot_res(nres_scwrl))=1
                  my_domain(my_tot_res(nres_scwrl))=iocc
                  my_bfactr(my_tot_res(nres_scwrl))=beta

C need to write the correct residue number - it'll be "i" (the residue
C number in the desired sequence) plus the offset determined earlier for
C this chain - this should ensure no repeats and no res #'s >9999...
                
c v2.1            write(num_of_res,'(i4)')mod(nres_scwrl,10000)
                  write(num_of_res,'(i4)')cur_offset(num_frags)+i

                  seq_scwrl(nres_scwrl:nres_scwrl)=
     &              seq_of_fasta(iii)(i:i)
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'D') 
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'E') 
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'H') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+0.5
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'K') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'R') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                endif
                if(strong(13:16).eq.' N  '.or.
     &             strong(13:16).eq.' CA '.or.
     &             strong(13:16).eq.' C  '.or.
     &             strong(13:16).eq.' O  ') then
                  strong(23:26)=num_of_res(1:4)
                  write(23,21)strong
                endif

C do the following if the matching res in the pdb is *similar*:

              elseif(i_in_pdb_fasta(i).eq.2.and.
     &               nres2.eq.my_nres2(i)) then
                if(strong(13:16).eq.' N  ') then
                  nres_scwrl=nres_scwrl+1
                  my_scwrl_res(num_tot_res_old+i)=nres_scwrl
                  my_tot_res(nres_scwrl)=num_tot_res_old+i

                  i_am_original(my_tot_res(nres_scwrl))=2
                  my_domain(my_tot_res(nres_scwrl))=iocc
                  my_bfactr(my_tot_res(nres_scwrl))=beta

C need to write the correct residue number - it'll be "i" (the residue
C number in the desired sequence) plus the offset determined earlier for
C this chain - this should ensure no repeats and no res #'s >9999...
                
c v2.1            write(num_of_res,'(i4)')mod(nres_scwrl,10000)
                  write(num_of_res,'(i4)')cur_offset(num_frags)+i

                  seq_scwrl(nres_scwrl:nres_scwrl)=
     &              seq_of_fasta(iii)(i:i)
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'D')
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'E')
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'H')
     &               qqq_tot_scwrl=qqq_tot_scwrl+0.5
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'K')
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'R')
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                endif
                if(strong(13:16).eq.' N  '.or.
     &             strong(13:16).eq.' CA '.or.
     &             strong(13:16).eq.' C  '.or.
     &             strong(13:16).eq.' O  ') then
                  strong(23:26)=num_of_res(1:4)
                  write(23,21)strong
                endif

C do the following if the matching res in the pdb is *identical*:

              elseif(i_in_pdb_fasta(i).eq.3.and.
     &               nres2.eq.my_nres2(i)) then
                if(strong(13:16).eq.' N  ') then
                  nres_scwrl=nres_scwrl+1
                  my_scwrl_res(num_tot_res_old+i)=nres_scwrl
                  my_tot_res(nres_scwrl)=num_tot_res_old+i

C note that if it's identical it likely won't get messed with in this
C run of ahemodel, so we need to not only store i_am_original (which we
C did above for non-identical matches), we also need to determine if it
C was model-built in a previous call to ahemodel - we can figure this
C out by asking if it's got zero in the beta factor field...
C we will use i_old_unstruc=1 later when we write out AHEMODEL...

                  i_am_original(my_tot_res(nres_scwrl))=3
                  my_domain(my_tot_res(nres_scwrl))=iocc
                  my_bfactr(my_tot_res(nres_scwrl))=beta

C need to write the correct residue number - it'll be "i" (the residue
C number in the desired sequence) plus the offset determined earlier for
C this chain - this should ensure no repeats and no res #'s >9999...
                
c v2.1            write(num_of_res,'(i4)')mod(nres_scwrl,10000)
                  write(num_of_res,'(i4)')cur_offset(num_frags)+i

                  n1=nres_scwrl ! just used to fit following lines into
C                               ! a single line of code...
                  if(seq_of_fasta(iii)(i:i).eq.'A') seq_scwrl(n1:n1)='a'
                  if(seq_of_fasta(iii)(i:i).eq.'C') seq_scwrl(n1:n1)='c'
                  if(seq_of_fasta(iii)(i:i).eq.'D') seq_scwrl(n1:n1)='d'
                  if(seq_of_fasta(iii)(i:i).eq.'E') seq_scwrl(n1:n1)='e'
                  if(seq_of_fasta(iii)(i:i).eq.'F') seq_scwrl(n1:n1)='f'
                  if(seq_of_fasta(iii)(i:i).eq.'G') seq_scwrl(n1:n1)='g'
                  if(seq_of_fasta(iii)(i:i).eq.'H') seq_scwrl(n1:n1)='h'
                  if(seq_of_fasta(iii)(i:i).eq.'I') seq_scwrl(n1:n1)='i'
                  if(seq_of_fasta(iii)(i:i).eq.'K') seq_scwrl(n1:n1)='k'
                  if(seq_of_fasta(iii)(i:i).eq.'L') seq_scwrl(n1:n1)='l'
                  if(seq_of_fasta(iii)(i:i).eq.'M') seq_scwrl(n1:n1)='m'
                  if(seq_of_fasta(iii)(i:i).eq.'N') seq_scwrl(n1:n1)='n'
                  if(seq_of_fasta(iii)(i:i).eq.'P') seq_scwrl(n1:n1)='p'
                  if(seq_of_fasta(iii)(i:i).eq.'Q') seq_scwrl(n1:n1)='q'
                  if(seq_of_fasta(iii)(i:i).eq.'R') seq_scwrl(n1:n1)='r'
                  if(seq_of_fasta(iii)(i:i).eq.'S') seq_scwrl(n1:n1)='s'
                  if(seq_of_fasta(iii)(i:i).eq.'T') seq_scwrl(n1:n1)='t'
                  if(seq_of_fasta(iii)(i:i).eq.'V') seq_scwrl(n1:n1)='v'
                  if(seq_of_fasta(iii)(i:i).eq.'W') seq_scwrl(n1:n1)='w'
                  if(seq_of_fasta(iii)(i:i).eq.'Y') seq_scwrl(n1:n1)='y' 
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'d') 
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'e') 
     &               qqq_tot_scwrl=qqq_tot_scwrl-1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'h') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+0.5
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'k') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                  if(seq_scwrl(nres_scwrl:nres_scwrl).eq.'r') 
     &               qqq_tot_scwrl=qqq_tot_scwrl+1.0
                endif
                strong(23:26)=num_of_res(1:4)
                write(23,21)strong
              endif
            enddo 
            goto 47
48          close(11)
            write(*,*)
            write(*,*)'# of pre-assigned domains in frag# ',num_frags
            write(*,*)' = ',num_domains(num_frags)
            write(*,*)
          enddo ! do the next chain in the pdb
        enddo   ! do the next pdb (there will only be one)
      enddo     ! do the next fasta

      close(25) ! close the final list of alignments...
      close(23) ! close the scwrl_inp.pdb also
      close(93) ! close the DSSP_vs_PSIPRED.txt file also...

C if use_SS_constraints and ifrst_time_thru=1 or 2 go back to 9119
C v1.4

      if(use_SS_constraints.and.ifrst_time_thru.le.2) then
        write(*,*)'about to go back and do alignments again'
        write(*,*)'ifrst_time_thru = ',ifrst_time_thru

C I think we need a rescue here if a pdb/chain has no match to an
C input sequence: we don't want its sequence length to become zero
C the following lines solve a problem with non-aligned fastas ending up
C being truncated to nothing and causing a problem at line 1951, i.e.
C at lines that look like:

C         write(12,my_fmt)seq_of_pdb_chn(nnn,mmm)
C    &                 (1:len_of_pdb_chn(nnn,mmm))

        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            inotdone=1
            do iii=1,num_seqs
              if(ialign(iii,nnn,mmm).eq.1) inotdone=0
            enddo
            if(inotdone.eq.1) 
     &        len_of_pdb_chn_fin(nnn,mmm)=len_of_pdb_chn(nnn,mmm)
          enddo
        enddo
 
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            write(*,9120)nnn,mmm,len_of_pdb_chn(nnn,mmm),
     &                           len_of_pdb_chn_fin(nnn,mmm)
9120        format('old & new template lengths ',4i8)
            if(truncate_template.and.ifrst_time_thru.eq.1) then
              if(len_of_pdb_chn(nnn,mmm).ne.
     &           len_of_pdb_chn_fin(nnn,mmm)) then
                write(*,*)'WARNING!! - we will truncate the template!'
                write(*,*)'this should be a good thing not a bad thing'
                write(*,*)'but be aware that we may not want to do it'
                write(*,*)
              endif
            elseif(ifrst_time_thru.eq.1) then
              if(len_of_pdb_chn(nnn,mmm).ne.
     &           len_of_pdb_chn_fin(nnn,mmm)) then
                write(*,*)'NOTE!!! - we could choose to truncate this'
                write(*,*)'template if we set truncate_template to yes'
                write(*,*)
              endif
            endif
          enddo
        enddo
ccc     stop
        goto 9119
      elseif(use_SS_constraints.and.ifrst_time_thru.eq.3) then
        write(*,*)'no need to go  back and do alignments again'
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            write(*,9120)nnn,mmm,len_of_pdb_chn(nnn,mmm),
     &                           len_of_pdb_chn_fin(nnn,mmm)
          enddo
        enddo
      endif
      if(truncate_template.and.ifrst_time_thru.eq.3) then
        write(*,*)
        write(*,*)'writing alignment scores for each chain'
        write(*,*)
        write(*,*)'score1 = standard sequence alignment; no SS cons'
        write(*,*)'score2 = sequence alignment w/ SS cons full-length'
        write(*,*)'score3 = sequence alignment w/ SS cons truncated'
        write(*,*)
        write(*,*)'if template was much too long then score3 should '
        write(*,*)'be considerably better than score2, but not as '
        write(*,*)'good as score1 (this should always be the best as'
        write(*,*)'no constraints are placed on the alignment...'
        write(*,*)
        do iii=1,num_seqs
          do nnn=1,num_pdbs
            do mmm=1,num_chns(nnn)
              if(ialign(iii,nnn,mmm).eq.0) cycle
              write(*,9121)iii,falign(iii,nnn,mmm)(1:6),
     &                        salign1(iii,nnn,mmm),
     &                        salign2(iii,nnn,mmm),
     &                        salign3(iii,nnn,mmm)
9121          format('Alignment scores ',i4,' ',a6,' ',3f10.3)
            enddo
          enddo
        enddo
        write(*,*)
      endif


C let's now make a rosetta.alignment file together with a combined fasta

      open(unit=82,file='rosetta.fasta',status='unknown')
      open(unit=83,file='rosetta.alignment',status='unknown')
      open(unit=84,file='temporary.file',status='unknown')
      write(82,'(">target")')
      write(83,'(">target")')
      write(84,'(">rosetta.pdb")')

C first we'll figure out how many chains there'll be at the end - this
C will differ depending on whether we use rosetta_symmetry

      if(rosetta_symmetry) then
        num_chns_final=num_seqs
      else
        num_chns_final=0 
        do iii=1,num_seqs
          do nnn=1,num_pdbs
            do mmm=1,num_chns(nnn)
              if(ialign(iii,nnn,mmm).eq.0) cycle
              num_chns_final=num_chns_final+1
            enddo
          enddo
        enddo
      endif

      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)

            if(ialign(iii,nnn,mmm).eq.0) cycle

C read the alignment associated with this sequence

            flen=len(trim(falign(iii,nnn,mmm)))
            falign2(iii,nnn,mmm)(1:flen)=falign(iii,nnn,mmm)
            falign2(iii,nnn,mmm)(flen+1:flen+10)='.rewritten'
            flen=flen+10
            falign2(iii,nnn,mmm)(flen+1:flen+8)='.rosetta'
            flen=flen+8
            write(*,*)'check falign2 ',falign2(iii,nnn,mmm)(1:flen)
            iseqs_read=0
            open(unit=93,file=falign2(iii,nnn,mmm)(1:flen),
     &           status='unknown')
234         read(93,'(a100)',end=345)string100
            if(string100(1:1).eq.'>') then
              iseqs_read=iseqs_read+1
              goto 234
            endif
            if(iseqs_read.eq.1) then
              write(83,'(a100)')string100
            elseif(iseqs_read.eq.2) then
              write(84,'(a100)')string100
            endif 
            goto 234
345         close(93)

C read the fasta associated with this sequence 

            fasta_tmp='XXXXXX.fasta'
            fasta_tmp(1:6)=falign2(iii,nnn,mmm)(1:6)
            open(unit=93,file=fasta_tmp,status='unknown')
235         read(93,'(a60)',end=346)string100(1:60)
            if(string100(1:1).eq.'>') goto 235
            write(82,'(a60)')string100(1:60)
            goto 235
346         close(93)

C if we're using rosetta_symmetry, then now that we've found a match 
C for this sequence we can skip to the next sequence - we only need 
C one copy in rosetta.alignment and rosetta.fasta when using symmetry
C write on a slash to rosetta.fasta before returning

            if(rosetta_symmetry) then
              if(iii.lt.num_seqs) write(82,'("/")')
              goto 123
            endif

C otherwise, if not using rosetta_symmetry then we still need to
C write on a slash to rosetta.fasta if we still have more fastas to add

            num_chns_tmp=num_chns_tmp+1
            if(num_chns_tmp.lt.num_chns_final) write(82,'("/")')

          enddo
        enddo
123     continue
      enddo
      close(82)
      close(83)
      close(84)

C make sure to add on the template's information into the final
C rosetta.alignment file

      call system('cat temporary.file >> rosetta.alignment')

C now write out TEMPLATE_rosetta.pdb in the chain-reordered way - note
C that we write out *all* chains regardless of symmetry being used

      open(unit=85,file='TEMPLATE_rosetta.pdb',status='unknown')
      num_chns_final=0
      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)

            if(ialign(iii,nnn,mmm).eq.0) cycle

C update total number of chains to be written to TEMPLATE_rosetta.pdb

            num_chns_final=num_chns_final+1

C read the input pdb and write out this chain to TEMPLATE_rosetta.pdb
C note that this final TEMPLATE_rosetta.pdb doesn't need SEQRES entries
C or anything like that...

            open(unit=93,file=fname_pdb_chn(nnn,mmm),status='unknown')
347         read(93,'(a66)',end=348)string66
            if(string66(1:4).ne.'ATOM') goto 347
            ijk=mod(num_chns_final,chainIDmax)
            if(ijk.eq.0) ijk=chainIDmax
            string66(22:22)=blah(ijk:ijk)
            write(85,'(a66)')string66
            goto 347
348         close(93)

          enddo
        enddo
      enddo

C add on any HETATM entries to the TEMPLATE_rosetta.pdb also
C v3.7 do the same for scwrl - put them in a scwrl_frame.pdb
C note that we use nscwrl_hets to determine whether to invoke later

C v4.1 it seems like scwrl must have some kind of in-built limit on the
C number of HETATM entries in scwrl_frame.pdb - it's continually
C crashing (very non-elegantly) when the #entries > 32500 or so - so
C we'll only add HETATM entries if they're within some cutoff of
C pre-existing ATOM entries - to determine this quickly we need to put
C the ATOM entries onto a grid...

      nat_tmp=0
      do nnn=1,num_pdbs
        open(unit=93,file=fname_pdb(nnn),status='unknown')
417     read(93,'(a66)',end=418)string66
        if(string66(1:4).eq.'ATOM') nat_tmp=nat_tmp+1
        goto 417
418     close(93)
      enddo

      allocate(xs(1:nat_tmp))
      allocate(ys(1:nat_tmp))
      allocate(zs(1:nat_tmp))

      nat_tmp=0
      do nnn=1,num_pdbs
        open(unit=93,file=fname_pdb(nnn),status='unknown')
517     read(93,'(a66)',end=518)string66
        if(string66(1:4).eq.'ATOM') then
          nat_tmp=nat_tmp+1
          read(string66,'(30x,3f8.3)')xt,yt,zt
          xs(nat_tmp)=xt
          ys(nat_tmp)=yt
          zs(nat_tmp)=zt
        endif
        goto 517
518     close(93)
      enddo

      xmin= 999999.9
      ymin= 999999.9
      zmin= 999999.9
      xmax=-999999.9
      ymax=-999999.9
      zmax=-999999.9
      do m=1,nat_tmp
        xmin=min(xmin,xs(m))
        ymin=min(ymin,ys(m))
        zmin=min(zmin,zs(m))
        xmax=max(xmax,xs(m))
        ymax=max(ymax,ys(m))
        zmax=max(zmax,zs(m))
      enddo
      xlen=xmax-xmin+10.0
      ylen=ymax-ymin+10.0
      zlen=zmax-zmin+10.0
      xinv=1.0/10.0 ! assume cell size of 10A
      yinv=1.0/10.0 
      zinv=1.0/10.0 
      nx=int(xlen/10.0)+1
      ny=int(ylen/10.0)+1
      nz=int(zlen/10.0)+1
      allocate(temp_grid(1:nx,1:ny,1:nz))
      temp_grid=0
      ninbox_max=0

      write(*,*)'xmin etc ',xmin,ymin,zmin
      write(*,*)'xmax etc ',xmax,ymax,zmax
      write(*,*)'xlen etc ',xlen,ylen,zlen
      write(*,*)'xinv etc ',xinv,yinv,zinv
      write(*,*)'nx   etc ',nx,ny,nz

      do m=1,nat_tmp
        xt=xs(m)
        yt=ys(m)
        zt=zs(m)
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              ninbox_max=max(ninbox_max,temp_grid(i1,j1,k1))
            enddo
          enddo
        enddo
      enddo
 
      write(*,*)'for HETATM check, ninbox_max = ',ninbox_max
      write(*,*)'allocating memory for grid'

      allocate(xatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))
      allocate(yatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))
      allocate(zatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))

      temp_grid=0

      do m=1,nat_tmp
        xt=xs(m)
        yt=ys(m)
        zt=zs(m)
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              xatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=xt
              yatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=yt
              zatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=zt
            enddo
          enddo
        enddo
      enddo

      write(*,*)'done putting everyone on the grid'

      nscwrl_hets_tot=0 ! includes *all* hetatm entries
      nscwrl_hets=0     ! only includes those within 10A of ATOMs
      open(unit=99,file='scwrl_frame.pdb',status='unknown')
      do nnn=1,num_pdbs
        open(unit=93,file=fname_pdb(nnn),status='unknown')
447     read(93,'(a66)',end=448)string66
        if(string66(1:4).ne.'HETA') goto 447
        if(string66(18:20).eq.'DUM') goto 447 ! skip OPM atoms
        nscwrl_hets_tot=nscwrl_hets_tot+1

C now determine if it's within 10A of an ATOM entry

        read(string66,'(30x,3f8.3)')xt,yt,zt

C assume it isn't for now (iwant=0)

        iwant=0
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1

C HETATMs might be off the grid (which only fits ATOM entries + 10A)
C so be prepared to skip this HETATM...

        if(i1.lt.1)  goto 447
        if(j1.lt.1)  goto 447
        if(k1.lt.1)  goto 447
        if(i1.gt.nx) goto 447
        if(j1.gt.ny) goto 447
        if(k1.gt.nz) goto 447

        do l1=1,temp_grid(i1,j1,k1)
          xu=xatm_grid(l1,i1,j1,k1) 
          yu=yatm_grid(l1,i1,j1,k1) 
          zu=zatm_grid(l1,i1,j1,k1) 
          dist2=(xt-xu)**2+
     &          (yt-yu)**2+
     &          (zt-zu)**2
          if(dist2.le.36.0) then ! set 6A as the limit
            iwant=1
            exit
          endif
        enddo

        if(iwant.eq.0) goto 447

C if we got here then we want it...

        write(85,'(a66)')string66

C now determine element for the scwrl_frame.pdb entry - default to
C assuming carbon...

        char1='C'
        if(string66(14:14).eq.'N') char1='N'
        if(string66(14:14).eq.'O') char1='O'
        if(string66(14:14).eq.'H') char1='H'
        if(string66(14:14).eq.'P') char1='P'
        if(string66(14:14).eq.'S') char1='S'
        if(string66(14:14).eq.'F') char1='F'
        write(99,'(a66,11x,a)')string66,char1
        nscwrl_hets=nscwrl_hets+1

        goto 447
448     close(93)
      enddo

      close(85)
      close(99)

      write(*,*)
      write(*,*)'#total HETATM entries = ',nscwrl_hets_tot
      write(*,*)'#within  6A of ATOMs  = ',nscwrl_hets
      write(*,*)

      deallocate(xs)
      deallocate(ys)
      deallocate(zs)
      deallocate(temp_grid)
      deallocate(xatm_grid)
      deallocate(yatm_grid)
      deallocate(zatm_grid)

C v2.7 regardless of align_mode we write out all statistics to fort.99

      if(align_mode.ge.1) then
        write(*,*)
        write(*,*)'now time to write out alignment stats to '
        write(*,*)'alignment_stats.txt'
        write(*,*)
        open(unit=99,file='alignment_stats.txt',status='unknown')
        write(99,'("REMARK")')
        do n=1,num_chns_tot ! v1.4
          write(99,843)n,chain_nam(n),nres_totl(n),nres_iden(n),
     &                  nres_posi(n)+nres_iden(n),
     &                  nres_alin(n),nres_gaps(n)
843       format('REMARK: chain # ',i3,' name: ',a6,' #res: ',i5,
     &          ' #iden: ',i6,' #iden+pos: ',i6,
     &          ' #other-align: ',i6,' #unalign: ',i6)
        enddo
        write(99,'("REMARK")')
        do n=1,num_chns_tot ! v1.4
          f1=100.0*real(nres_iden(n))/real(nres_totl(n))
          f2=100.0*real((nres_iden(n)+nres_posi(n)))/real(nres_totl(n))
          f3=100.0*real(nres_alin(n))/real(nres_totl(n))
          f4=100.0*real(nres_gaps(n))/real(nres_totl(n))
          write(99,844)n,chain_nam(n),nres_totl(n),f1,f2,f3,f4
844       format('REMARK: chain # ',i3,' name: ',a6,' #res: ',i5,
     &          ' %iden: ',f6.2,' %iden+pos: ',f6.2,
     &          ' %other-align: ',f6.2,' %unalign: ',f6.2)
        enddo
        write(99,'("REMARK")')
        do n=1,num_chns_tot ! v1.4

C note we here subtract off any preceding or succeeding residues to get
C the total number of unaligned template residues *within* alignment...
C note that we can only do this if they are not equal to -1

          if(nres_befr(n).gt.0) nres_gaps(n)=nres_gaps(n)-nres_befr(n)
          if(nres_aftr(n).gt.0) nres_gaps(n)=nres_gaps(n)-nres_aftr(n)
          
          write(99,845)n,chain_nam(n),nres_totl(n),nres_gaps(n),
     &                nres_befr(n),nres_aftr(n)
845       format('REMARK: chain # ',i3,' name: ',a6,' #res: ',i5,
     &          ' target-res #within-unalign ',i5,
     &          ' #prior ',i5,' #after ',i5)
        enddo
        write(99,'("REMARK")')
        write(*,*)'if #unalign is large you will be building large '
        write(*,*)'or many loops using loopy'
        write(*,*)
        write(*,*)'if #prior or #after are large you may need to use'
        write(*,*)'rosetta to add in secondary structure elements'
        write(*,*)
        do n=1,num_chns_tot ! v1.4

C note we here subtract off any preceding or succeeding residues to get
C the total number of unaligned template residues *within* alignment...

          if(nres_tmplte_befr(n).gt.0) nres_tmplte_gaps(n)=
     &       nres_tmplte_gaps(n)-nres_tmplte_befr(n)
          if(nres_tmplte_aftr(n).gt.0) nres_tmplte_gaps(n)=
     &       nres_tmplte_gaps(n)-nres_tmplte_aftr(n)

          write(99,846)n,chain_nam(n),nres_totl(n),nres_tmplte_gaps(n),
     &                        nres_tmplte_befr(n),nres_tmplte_aftr(n)
846       format('REMARK: chain # ',i3,' name: ',a6,' #res: ',i5,
     &          ' tmplte-res #within-unalign ',i5,
     &          ' #prior ',i5,' #after ',i5)
        enddo
        write(*,*)
        write(*,*)'if #unalign or #prior or #after are large you may'
        write(*,*)'get a better alignment by truncating the template'
        write(*,*)'C-terminal tails can be automatically removed using'
        write(*,*)'truncate_template = "yes"...'
        write(*,*)
        write(*,*)'finally, here is a list of all HETATM residue types'
        write(*,*)'found in the template - you may want these to be '
        write(*,*)'omitted from clash checking'
        write(*,*)
        write(99,'("REMARK")')
        close(99)
        if(num_hets_found.gt.0) then
          write(my_fmt3,'(a4,i0,a8)'),'(i6,',num_hets_found,'(1x,a3))'
          write(*,my_fmt3)num_hets_found,
     &                  (nam_hets_found(i2),i2=1,num_hets_found)
          write(*,*)
        endif
        if(align_mode.eq.2) then
          write(*,*)
          write(*,*)'all alignments made and since align_mode = 2 '
          write(*,*)'we can now quit'
          write(*,*)
          stop
        endif
      endif

C v4.1 following are old comments prior to figuring out that we can use a
C "frame" file with scrwl - we'll keep these comments in for now just in
C case we subsequently find problems with the frame approach...
C "-f scwrl_frame.pdb"      

CCC we tried the following but it didn't work:

CCC we want scwrl to see the HETATM entries when it places rotamers - we
CCC thought we could fake it out by adding these atoms as fake GLY
CCC residues (cycling through atom names N, CA, C O) but no matter how I
CCC tried to do this, scwrl always complained about the added 'GLY'
CCC residues - it seems that it knows that they're not real residues and
CCC complains as follows:

CCC Error: Backbone's peristency failure in residue "C 861 "

CCC I tried renumbering residues, changing chain name etc but nothing
CCC fixed this problem - I therefore think that we should just go ahead
CCC and put up with this for now

CCC so, instead of using code at this point to fake out scwrl, I'm going
CCC to take the same code and move it to the point at which we call loopy
CCC so that loopy will hopefully try to avoid HETATM entries...

CCC UPDATE - I think we should look again at the options for scwrl as I'm
CCC pretty sure that it's supposed to be able to handle things like this
CCC so:

      write(my_fmt,'(a,i0,a)'),'(a',nres_scwrl,')'
      write(24,my_fmt)seq_scwrl(1:nres_scwrl)

      open(unit=25,file='seq_final.fasta',status='unknown')
      write(my_fmt,'(a,i0,a)'),'(a',num_tot_res,')'
      write(25,my_fmt)seq_tot_res(1:num_tot_res)

      close(24)
      close(25) 

C v2.0 we are now using an indirect method of counting residues in the
C scwrl_out.pdb in an attempt to make it work for >10,000 residues - we
C should therefore put in a check that the nres_scwrl that we know at
C this point matches with the one that we figure out later...

      nres_scwrl_inp=nres_scwrl

C 2023 reset all loops to zero if skip_loopy

      if(skip_loopy) num_tot_loops=0

C at this point we can read in any overrides for N- and C-terminal tails
C read "N" or "C", then the chain#, then the new (chain-local) res# for 
C the start of the tail to rebuild - note for N-terminal tails we only
C need to update iend and for C-terminal tails we only need to update
C ibeg - since we only input the chain-local res# we need to also figure
C out the offset (ndiff) needed to get the global res#

      irebuild_nters=0 ! set to 1 if this Nter tail is rebuilt
      irebuild_cters=0

      num_tot_nters_before_rebuild=num_tot_nters
      num_tot_cters_before_rebuild=num_tot_cters

      write(*,*)'initial #NTERS = ',num_tot_nters
      write(*,*)'initial #CTERS = ',num_tot_cters

      if(rebuild_tails) then
        write(*,*)
        write(*,*)'we wish to rebuild/redefine N- and C-terminal tails'
        write(*,*)'so now read "tail_rebuilds.txt"'
        write(*,*)
        write(*,*)'format is (1) "N" or "C", (2) # of chain '
        write(*,*)'(3) # of res where now to start tail'
        write(*,*)
        open(unit=25,file='tail_rebuilds.txt',status='unknown')
1101    read(25,*,end=1102)atail,itail,jtail

C we look to see if the tail to be rebuilt was already designated as a
C tail (if so, ifound=1) - need to also handle case where it's a new one

C note that here we make use of:

C my_frst_res(1:# total chains) - frst res# in total list
C my_last_res(1:# total chains) - last res# in total list
C my_totl_res(1:# total chains) - # of res in this particular chain

        if(itail.gt.num_frags) then
          write(*,*)
          write(*,*)'tail_rebuilds.txt list a chain # that '
          write(*,*)'exceeds the total number of chains (num_frags) '
          write(*,*)'in the pdb file ',num_frags
          write(*,*)'this is an error so quitting :('
          write(*,*)'please edit tail_rebuilds.txt to fix this issue'
          write(*,*) 
          stop
        endif

        if(atail.eq.'N') then
          ifound=0
          do n=1,num_tot_nters
            if(itail.eq.ifrg_tot_nters(n)) then
              ndiff=iend_tot_nters(n)-iend_loc_nters(n)
              iend_loc_nters(n)=jtail
              iend_tot_nters(n)=jtail+ndiff
              irebuild_nters(n)=1
              ifound=1
            endif
          enddo
          if(ifound.eq.0) then
            num_tot_nters=num_tot_nters+1
            irebuild_nters(num_tot_nters)=1
            ibeg_loc_nters(num_tot_nters)=1
            iend_loc_nters(num_tot_nters)=jtail
            ifrg_tot_nters(num_tot_nters)=itail
            ichn_tot_nters(num_tot_nters)=cur_chainID(itail)
            ibeg_tot_nters(num_tot_nters)=my_frst_res(itail)
            iend_tot_nters(num_tot_nters)=my_frst_res(itail)+jtail-1
          endif
        elseif(atail.eq.'C') then
          ifound=0
          do n=1,num_tot_cters
            if(itail.eq.ifrg_tot_cters(n)) then
              ndiff=ibeg_tot_cters(n)-ibeg_loc_cters(n)
              ibeg_loc_cters(n)=jtail
              ibeg_tot_cters(n)=jtail+ndiff
              irebuild_cters(n)=1
              ifound=1
            endif
          enddo
          if(ifound.eq.0) then
            num_tot_cters=num_tot_cters+1
            irebuild_cters(num_tot_cters)=1
            ibeg_loc_cters(num_tot_cters)=jtail
            iend_loc_cters(num_tot_cters)=my_totl_res(itail)
            ifrg_tot_cters(num_tot_cters)=itail
            ichn_tot_cters(num_tot_cters)=cur_chainID(itail)
            ibeg_tot_cters(num_tot_cters)=my_frst_res(itail)+jtail-1
            iend_tot_cters(num_tot_cters)=my_last_res(itail)
          endif
        endif 
        goto 1101
1102    close(25)

      elseif(rebuild_auto) then
        write(*,*)
        write(*,*)'we wish to rebuild/redefine N- and C-terminal tails'
        write(*,*)'so now read "auto_rebuilds.txt"'
        write(*,*)
        write(*,*)'format is: '
        write(*,*)'(1) # of chain '
        write(*,*)'(2) domain # of Nterminal stationary part'
        write(*,*)'(3) domain # of Cterminal stationary part'
        write(*,*)
        open(unit=25,file='auto_rebuilds.txt',status='unknown')
1103    read(25,*,end=1104)itail,idomN,idomC

        if(itail.gt.num_frags) then
          write(*,*)
          write(*,*)'auto_rebuilds.txt list a chain # that '
          write(*,*)'exceeds the total number of chains (num_frags) '
          write(*,*)'in the pdb file ',num_frags
          write(*,*)'this is an error so quitting :('
          write(*,*)'please edit tail_rebuilds.txt to fix this issue'
          write(*,*) 
          stop
        endif

C first, we write out (for checking purposes only) the domains found in
C this chain, together with their first residue number

        itmp=0
        do j=my_frst_res(itail),my_last_res(itail)
          if(i_in_dom_fasta_tot(j).gt.itmp) then
            write(*,*)'SANITY CHECK for chain ',itail,
     &                ' we have a new domain ',i_in_dom_fasta_tot(j),
     &                ' starting at (total) res# ',j
            itmp=i_in_dom_fasta_tot(j)
          elseif(i_in_dom_fasta_tot(j).lt.itmp) then
            write(*,*)'SANITY CHECK for chain ',itail,
     &                ' that previous domain              ',
     &                ' now ends at (total) res# ',j
            itmp=i_in_dom_fasta_tot(j)
          endif
        enddo

C now find the first residue of "idomN"...
         
        do j=my_frst_res(itail),my_last_res(itail)
          if(i_in_dom_fasta_tot(j).eq.idomN) then
            jtailN=j-my_frst_res(itail)+1
            write(*,*)'jtailN initial value ',jtailN

C now continue to count forwards till we find an unstructured residue -
C this way, we keep loop atoms that might have been added in a
C previous iteration and don't try to rebuild them - note that the first
C unstructured residue that we find is the one that we want as the last
C residue of the tail to be built

            ifound=0
            do k=j,1,-1
              if(i_in_dom_fasta_tot(k).ne.idomN.and.
     &           i_in_dom_fasta_tot(k).ne.0) then
                jtailN=k-my_frst_res(itail)+1
                ifound=1
                exit
              endif
            enddo
            if(ifound.eq.0) jtailN=0 ! I think we have to do this
            write(*,*)'jtailN final   value ',jtailN

            exit
          endif
       
        enddo

C ...and find the last residue of "idomC"

        do j=my_last_res(itail),my_frst_res(itail),-1
          if(i_in_dom_fasta_tot(j).eq.idomC) then
            jtailC=j-my_frst_res(itail)+1
            write(*,*)'jtailC initial value ',jtailC

C now continue to count backwards till we find an unstructured residue -
C this way, we keep loop atoms that might have been added in a
C previous iteration and don't try to rebuild them - note that the first
C unstructured residue that we find is the one that we want as the last
C residue of the tail to be built

            ifound=0
            do k=j,my_last_res(itail)
              if(i_in_dom_fasta_tot(k).ne.idomC.and.
     &           i_in_dom_fasta_tot(k).ne.0) then
                jtailC=k-my_frst_res(itail)+1
                ifound=1
                exit
              endif
            enddo
            if(ifound.eq.0) jtailC=my_totl_res(itail)+1 ! necessary
            write(*,*)'jtailC final   value ',jtailC

            exit
          endif
        enddo

C and use these to determine if we need to add a new tail or not - and
C where the tails should stop - note that the values thus determined are
C preliminary ones that will probably be overwritten later - but we do
C this now for convenience - we can reuse the code used above for
C "rebuild_tails" - to find the final values determined for tail termini 
C just look for "rebuild_auto" further down the code...

C note that what we are most interested in defining at this point in the
C code are: iend_loc_nters & ibeg_loc_cters (and their tot equivalents)
C later on in the end code we are most interested in defining instead
C the following: ibeg_loc_nters & iend_loc_cters...

ccc     write(*,*)'jtailN & jtailC = ',jtailN,jtailC

C from here the code is identical to that used in rebuild_tails...

C we look to see if the tail to be rebuilt was already designated as a
C tail (if so, ifound=1) - need to also handle case where it's a new one

C note that here we make use of:

C my_frst_res(1:# total chains) - frst res# in total list
C my_last_res(1:# total chains) - last res# in total list
C my_totl_res(1:# total chains) - # of res in this particular chain

C if we decided earlier that this tail was fully built we shouldn't
C attempt to do anything here...

        if(jtailN.eq.0) goto 6006
        ifound=0
        do n=1,num_tot_nters
          if(itail.eq.ifrg_tot_nters(n)) then
            ndiff=iend_tot_nters(n)-iend_loc_nters(n)
            iend_loc_nters(n)=jtailN
            iend_tot_nters(n)=jtailN+ndiff
            irebuild_nters(n)=1
            ifound=1
          endif
        enddo
        if(ifound.eq.0) then
          num_tot_nters=num_tot_nters+1
          irebuild_nters(num_tot_nters)=1
          ibeg_loc_nters(num_tot_nters)=1
          iend_loc_nters(num_tot_nters)=jtailN
          ifrg_tot_nters(num_tot_nters)=itail
          ichn_tot_nters(num_tot_nters)=cur_chainID(itail)
          ibeg_tot_nters(num_tot_nters)=my_frst_res(itail)
          iend_tot_nters(num_tot_nters)=my_frst_res(itail)+jtailN-1
        endif
6006    continue

C if we decided earlier that this tail was fully built we shouldn't
C attempt to do anything here...

        if(jtailC.eq.my_totl_res(itail)+1) goto 7007
        ifound=0
        do n=1,num_tot_cters
          if(itail.eq.ifrg_tot_cters(n)) then
            ndiff=ibeg_tot_cters(n)-ibeg_loc_cters(n)
            ibeg_loc_cters(n)=jtailC
            ibeg_tot_cters(n)=jtailC+ndiff
            irebuild_cters(n)=1
            ifound=1
          endif
        enddo
        if(ifound.eq.0) then
          num_tot_cters=num_tot_cters+1
          irebuild_cters(num_tot_cters)=1
          ibeg_loc_cters(num_tot_cters)=jtailC
          iend_loc_cters(num_tot_cters)=my_totl_res(itail)
          ifrg_tot_cters(num_tot_cters)=itail
          ichn_tot_cters(num_tot_cters)=cur_chainID(itail)
          ibeg_tot_cters(num_tot_cters)=my_frst_res(itail)+jtailC-1
          iend_tot_cters(num_tot_cters)=my_last_res(itail)
        endif
7007    continue

C now go back and read another line of auto_rebuilds.txt

        goto 1103
1104    close(25)

      endif ! close the if rebuild_tails/rebuild_auto

C note that we're not *quite* done with analyzing loops: there is a
C chance that some loops might have become Nter or Cter tails when they
C were expanded to contain min_loop_size residues - if so, we need to
C make sure they get put in the tail list and 'removed' from loop list -
C note that since removing is a pain in the ass we'll simply set
C ibeg_loc_loops (and iend_loc_loops) to -999 so that we can easily
C identify them and skip them when we reach the loopy part...
C note that below we also make use of ifrg_tot_loops, ifrg_tot_nters etc
C - these are the chain # of this particular chain

C note that we use i_done_res to ensure that structured parts of the
C "tail" are passed to loopy *unless* they are in to-be-built loops

      do i=1,num_tot_nters
        write(*,*)'NTERS ',i,ibeg_tot_nters(i),iend_tot_nters(i)
        do j=ibeg_tot_nters(i),iend_tot_nters(i)

C v2.0 - following lines now commented out - seems unnecessary (and
C wrong) to set i_am_original=0 for all residues in prior tails

ccc       if(i.le.num_tot_nters_before_rebuild) then
ccc         i_am_original(j)=0 ! used to ensure we build it later
ccc       endif

          i_done_res(j)=1    ! ensures res passed from scwrl to loopy
C                            ! needed when we use rebuild tails option
          do k=1,num_tot_loops
            if(j.ge.ibeg_tot_loops(k).and.j.le.iend_tot_loops(k)) then
              i_done_res(j)=0 ! change mind if this res is in a loop
            endif
          enddo
        enddo
      enddo

      do i=1,num_tot_cters
        write(*,*)'CTERS ',i,ibeg_tot_cters(i),iend_tot_cters(i)
        do j=ibeg_tot_cters(i),iend_tot_cters(i)

C v2.0 - following lines now commented out - seems unnecessary (and
C wrong) to set i_am_original=0 for all residues in prior tails

ccc       if(i.le.num_tot_cters_before_rebuild) then
ccc         i_am_original(j)=0 ! used to ensure we build it later
ccc       endif

          i_done_res(j)=1    ! ensures res passed from scwrl to loopy
C                            ! needed when we use rebuild tails option
          do k=1,num_tot_loops
            if(j.ge.ibeg_tot_loops(k).and.j.le.iend_tot_loops(k)) then
              i_done_res(j)=0 ! change mind if this res is in a loop
            endif
          enddo
        enddo
      enddo

      do i=1,num_tot_loops

C make loops that are at beg of chain become new Nter tails - I think
C this is unlikely to happen but let's plan for it anyway

C UPDATE - I think we should only do this if there is not already a Nter
C tail, so check for that first...

        do j=1,num_tot_nters
          if(ifrg_tot_nters(j).eq.ifrg_tot_loops(i)) goto 456
        enddo

        if(ibeg_loc_loops(i).eq.1) then
          write(*,*)'loop #',i,' is a new Nter tail ',num_tot_nters
          num_tot_nters=num_tot_nters+1
          ifrg_tot_nters(num_tot_nters)=ifrg_tot_loops(i)
          ichn_tot_nters(num_tot_nters)=cur_chainID(ifrg_tot_loops(i))
          ibeg_tot_nters(num_tot_nters)=ibeg_tot_loops(i)
          iend_tot_nters(num_tot_nters)=iend_tot_loops(i)
          ibeg_loc_nters(num_tot_nters)=ibeg_loc_loops(i)
          iend_loc_nters(num_tot_nters)=iend_loc_loops(i)
          ibeg_loc_loops(i)=-999 ! exploit this later
          iend_loc_loops(i)=-999 ! exploit this later
          do j=ibeg_tot_nters(num_tot_nters),
     &         iend_tot_nters(num_tot_nters)
            i_am_original(j)=0 ! make sure set to zero
          enddo
        elseif(ibeg_loc_loops(i).lt.1.and.
     &         ibeg_loc_loops(i).ne.-999) then
          write(*,*)'error for beg loop i ',i,ibeg_loc_loops(i)
          stop
        endif

456     continue

C make loops that are at end of chain become new Cter tails - I think
C this is also unlikely to happen but let's plan for it anyway

C UPDATE - I think we should only do this if there is not already a Nter
C tail, so check for that first...

        do j=1,num_tot_cters
          if(ifrg_tot_cters(j).eq.ifrg_tot_loops(i)) goto 567
        enddo

        if(iend_loc_loops(i).eq.num_res_in_frag(ifrg_tot_loops(i))) then
          num_tot_cters=num_tot_cters+1
          write(*,*)'loop #',i,' is a new Cter tail ',num_tot_cters
          ifrg_tot_cters(num_tot_cters)=ifrg_tot_loops(i)
          ichn_tot_cters(num_tot_cters)=cur_chainID(ifrg_tot_loops(i))
          ibeg_tot_cters(num_tot_cters)=ibeg_tot_loops(i)
          iend_tot_cters(num_tot_cters)=iend_tot_loops(i)
          ibeg_loc_cters(num_tot_cters)=ibeg_loc_loops(i)
          iend_loc_cters(num_tot_cters)=iend_loc_loops(i)
          ibeg_loc_loops(i)=-999 ! exploit this later
          iend_loc_loops(i)=-999 ! exploit this later
          do j=ibeg_tot_cters(num_tot_cters),
     &         iend_tot_cters(num_tot_cters)
            i_am_original(j)=0 ! make sure set to zero
          enddo
        elseif(iend_loc_loops(i).gt.num_res_in_frag(ifrg_tot_loops(i)))
     &  then
          write(*,*)'error for end loop i ',i,iend_loc_loops(i)
          stop
        endif

567     continue

C make loops that encroach on existing Nter tails part of those tails...
C note that if the tail is not being rebuilt then we do normal code but
C if it is being rebuilt then we only "remove" the loop if it contains
C the final residue of the new tail - this allows us to retain other
C loops that might be built in a structured domain at the N-terminal
C side of the new tail position - note also that if we do "remove" the
C loop then we need to update the start position of the tail so that it
C builds in all unresolved residues - this is actually a neat feature
C because it means we don't need to provide an exactly-correct res# for
C the position of the rebuilt N-terminal tail: so long as we get it at
C least within the inter-domain loop the code should automatically
C figure out the correct starting point

        do j=1,num_tot_nters
          if(ifrg_tot_loops(i).ne.ifrg_tot_nters(j)) cycle
          if(irebuild_nters(j).eq.0) then
            if(ibeg_loc_loops(i).le.iend_loc_nters(j).and.
     &         ibeg_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' merges with Nter tail ',j

C if the loop ends further back than the Nter then update end of Nter

              if(iend_loc_nters(j).lt.iend_loc_loops(i)) then
                iend_tot_nters(j)=iend_tot_loops(i)
                iend_loc_nters(j)=iend_loc_loops(i)
              endif
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          elseif(irebuild_nters(j).eq.1) then
            if(ibeg_loc_loops(i).le.iend_loc_nters(j).and.
     &         iend_loc_loops(i).ge.iend_loc_nters(j).and.
     &         ibeg_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' merges with rebuilt Nter tail ',j
              write(*,*)'Nter tail no longer starts at local res # ',
     &                   iend_loc_nters(j)
              write(*,*)'Nter tail now starts at local res # ',
     &                   iend_loc_loops(i)
              iend_tot_nters(j)=iend_tot_loops(i)
              iend_loc_nters(j)=iend_loc_loops(i)
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          endif
        enddo

C make loops that encroach on existing Cter tails part of those tails...

        do j=1,num_tot_cters
          if(ifrg_tot_loops(i).ne.ifrg_tot_cters(j)) cycle
          if(irebuild_cters(j).eq.0) then
            if(iend_loc_loops(i).ge.ibeg_loc_cters(j).and.
     &         iend_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' merges with Cter tail ',j

C if the loop begs further fwrd than the Cter then update beg of Cter

              if(ibeg_loc_cters(j).gt.ibeg_loc_loops(i)) then
                ibeg_tot_cters(j)=ibeg_tot_loops(i)
                ibeg_loc_cters(j)=ibeg_loc_loops(i)
              endif
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          elseif(irebuild_cters(j).eq.1) then
            if(ibeg_loc_loops(i).le.ibeg_loc_cters(j).and.
     &         iend_loc_loops(i).ge.ibeg_loc_cters(j).and.
     &         iend_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' merges with rebuilt Cter tail ',j
              write(*,*)'Cter tail no longer starts at local res # ',
     &                   ibeg_loc_cters(j)
              write(*,*)'Cter tail now starts at local res # ',
     &                   ibeg_loc_loops(i)
              ibeg_tot_cters(j)=ibeg_tot_loops(i)
              ibeg_loc_cters(j)=ibeg_loc_loops(i)
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          endif
        enddo
      enddo

C now write out the final final loops and tails...
C *and* set i_am_original to loops to zero...

      do i=1,num_tot_loops
        if(ibeg_loc_loops(i).gt.0) then
          do j=ibeg_tot_loops(i),iend_tot_loops(i)
            i_am_original(j)=0 ! set all res in real loops to zero
          enddo

C v1.9 - for writing out the seq of each loop we'll use itl - this will
C be the ID of the sequence that this loop's chain matched to - we
C figure out itl in two stages: (1) first, set itl to the fragment # of
C this loop (i.e. its chain #), then (2) set itl to the sequence #
C matched to this fragment earlier on...

          itl=ifrg_tot_loops(i) ! just use for writing loop seq here
          itl=my_frags_seq_num(itl)

          write(*,778) i,ibeg_tot_loops(i),iend_tot_loops(i),
     &                   ichn_tot_loops(i),ifrg_tot_loops(i),
     &                   ibeg_loc_loops(i),iend_loc_loops(i),
     & seq_of_fasta(itl)(ibeg_loc_loops(i):ibeg_loc_loops(i)+2),
     & seq_of_fasta(itl)(iend_loc_loops(i)-2:iend_loc_loops(i))
778       format('FINAL FINAL LOOPS ',3i8,1x,a,1x,3i8,1x,
     &            a3,'...',a3)
        endif
      enddo

      do i=1,num_tot_nters
        write(*,978)i,ibeg_tot_nters(i),iend_tot_nters(i)
978     format('FINAL FINAL NTERS (prior to possible truncation) ',3i8)
      enddo

      do i=1,num_tot_cters
        write(*,979)i,ibeg_tot_cters(i),iend_tot_cters(i)
979     format('FINAL FINAL CTERS (prior to possible truncation) ',3i8)
      enddo

      if(truncate_tails) then
        write(*,*)
        write(*,*)'we wish to truncate N- or C-terminal tails'
        write(*,*)'so now read "tail_truncations.txt"'
        write(*,*)
        write(*,*)'format is (1) "N" or "C", (2) # of chain '
        write(*,*)'(3) # of chain-local res where now to truncate tail'
        write(*,*)
        open(unit=25,file='tail_truncations.txt',status='unknown')
        nfound=0
1201    read(25,*,end=1202)atail,itail,jtail
        if(atail.eq.'N') then
          ifound=0
          do n=1,num_tot_nters
            if(itail.eq.ifrg_tot_nters(n)) then
              ndiff=jtail-ibeg_loc_nters(n)
              ibeg_loc_nters(n)=jtail
              ibeg_tot_nters(n)=ibeg_tot_nters(n)+ndiff
              ifound=1
              nfound=nfound+1
            endif
          enddo
          if(ifound.eq.0) then
            write(*,*)'you are truncating a N-ter tail (#',itail,
     &                ') that does not appear to exist - this is fatal'
            write(*,*)'quitting now :('
            stop
          endif
        elseif(atail.eq.'C') then
          ifound=0
          do n=1,num_tot_cters
            if(itail.eq.ifrg_tot_cters(n)) then
              ndiff=iend_loc_cters(n)-jtail
              iend_loc_cters(n)=jtail
              iend_tot_cters(n)=iend_tot_cters(n)-ndiff
              ifound=1
              nfound=nfound+1
            endif
          enddo
          if(ifound.eq.0) then
            write(*,*)'you are truncating a C-ter tail (#',itail,
     &                ') that does not appear to exist - this is fatal'
            write(*,*)'quitting now :('
            stop
          endif
        endif
        goto 1201 ! read another tail
1202    close(25)
        if(nfound.eq.0) then
          write(*,*)'you set truncate_tails = "yes" but there was '
          write(*,*)'nothing in "tail_truncations.txt"'
          write(*,*)'quitting now :('
          stop
        endif

C we also may truncate tails if rebuild_auto (these are mutual exclus)
C here we have a crude recompilation option - if we uncomment the
C "dogshit" line then the code will skip the determination of new
C truncation points and will instead attempt to build the entire
C molecule using the rebuild_tails type of approach - this may work in
C relatively simple cases - if, however, we use the "rebuild_auto" if
C statement then the code will read "auto_rebuilds.txt" with "itail"
C being the identity of the chain, and "idomN" and "idomC" being the
C identities of the N-terminal-most and C-terminal-most *fixed* domains
C in the chain - once these have been read the code will determine where
C the tail preceding the domain preceding "idomN" ends, and where the 
C tail succeeding the domain succeeding "idomC" ends - it will then
C attempt to build to those points subsequently - the reasoning is that
C we should probably only attempt to rebuild one domain at a time in
C each of the N- and C-directions, then save the coordinates and rerun
C ahemodel iteratively until all desired domains have been added - we
C try to build the linker/tails preceding the new N-domain and
C succeeding the new C-domain since this gives us leeway in running the
C next iteration: if we can't add on a further domain then the code will
C gradually "eat away" at a previously-built loop rather than eating
C away directly into a structured domain...

      elseif(idogshit.eq.999) then
ccc   elseif(rebuild_auto) then

        open(unit=25,file='auto_rebuilds.txt',status='unknown')
1203    read(25,*,end=1204)itail,idomN,idomC

C find how many domains there are in total in this chain

        maxdom=0
        do j=my_frst_res(itail),my_last_res(itail)
          maxdom=max(maxdom,i_in_dom_fasta_tot(j))
        enddo

C find the first "total" residue of domain idomN, jtailN
CCC find the first chain-local residue of domain idomN, jtailN

        do j=my_frst_res(itail),my_last_res(itail)
          if(i_in_dom_fasta_tot(j).eq.idomN) then
            jtailN=j
ccc         jtailN=j-my_frst_res(itail)+1
            exit
          endif
        enddo

C find the last "total" residue of domain idomC, jtailC 
CCC find the last chain-local residue of domain idomC, jtailC

        do j=my_last_res(itail),my_frst_res(itail),-1
          if(i_in_dom_fasta_tot(j).eq.idomC) then
            jtailC=j
ccc         jtailC=j-my_frst_res(itail)+1
            exit
          endif
        enddo

C if there are at least two domains to build prior to idomN then count
C backwards from jtailN (defined above) till we find the last residue of
C domain idomN-2, then add one so that we have the last loop residue
C note that jtailN_final is a chain-local number whereas jtailN is not

        if(idomN.gt.2) then

          do j=jtailN,1,-1
            if(i_in_dom_fasta_tot(j).eq.idomN-2) then
              jtailN_final=j-my_frst_res(itail)+1
              jtailN_final=jtailN_final+1 ! want res up to the domain
              exit
            endif
          enddo

C then use this number to reset ibeg_loc_nters and ibeg_tot_nters...

          do n=1,num_tot_nters
            if(itail.eq.ifrg_tot_nters(n)) then
              ndiff=jtailN_final-ibeg_loc_nters(n)
              ibeg_loc_nters(n)=jtailN_final
              ibeg_tot_nters(n)=ibeg_tot_nters(n)+ndiff
            endif
          enddo

        endif

C if there are at least two domains to build after idomC then count
C forwards from jtailC (defined above) till we find the first residue of
C domain idomC+2, then subtract one so that we have the last loop residue
C note that jtailC_final is a chain-local number whereas jtailC is not

        if(idomC.lt.maxdom-1) then

          do j=jtailC,my_last_res(itail)
            if(i_in_dom_fasta_tot(j).eq.idomC+2) then
              jtailC_final=j-my_frst_res(itail)+1
              jtailC_final=jtailC_final-1 ! want res up to the domain
              exit
            endif
          enddo

C then use this number to reset iend_loc_cters and iend_tot_cters...

          do n=1,num_tot_cters
            if(itail.eq.ifrg_tot_cters(n)) then
              ndiff=iend_loc_cters(n)-jtailC_final
              iend_loc_cters(n)=jtailC_final
              iend_tot_cters(n)=iend_tot_cters(n)-ndiff
            endif
          enddo

        endif

        goto 1203
1204    close(25)

      endif ! close the if for tail truncation

C finally write out the new tails 
 
      do i=1,num_tot_nters
c       write(*,878)i,ibeg_tot_nters(i),iend_tot_nters(i)
c878    format('FINAL FINAL NTERS (after possible truncation) ',3i8)
        write(*,878)i,ifrg_tot_nters(i),
     &                ibeg_tot_nters(i),iend_tot_nters(i)
878     format('FINAL FINAL NTERS (after possible truncation) '
     &         ' nter# ',i4,' chn# ',i4,' beg ',i6,' end ',i6)
      enddo

      do i=1,num_tot_cters
c       write(*,879)i,ibeg_tot_cters(i),iend_tot_cters(i)
c879    format('FINAL FINAL CTERS (after possible truncation) ',3i8)
        write(*,879)i,ifrg_tot_cters(i),
     &                ibeg_tot_cters(i),iend_tot_cters(i)
879     format('FINAL FINAL CTERS (after possible truncation) '
     &         ' cter# ',i4,' chn# ',i4,' beg ',i6,' end ',i6)
      enddo

      do i=1,num_tot_res
        write(*,919)i,my_scwrl_res(i),res_tot_res(i),
     &              i_am_original(i),i_in_dom_fasta_tot(i)
919     format('Check Scwrl Residues ',5i8)
      enddo

C v2.3 also check that if we have the chainID and the resID from a
C loopy_inp.pdb we can convert back to a total res#...

      do nyyy=1,num_frags
        do nj=1,chainIDmax
          if(cur_chainID(nyyy).eq.blah(nj:nj)) then
            nxxx=nj
            exit
          endif
        enddo
        do nz=1,num_res_in_frag(nyyy)
          write(*,*)'CHECK cur_chainID ',cur_chainID(nyyy),
     &              ' res# ',nz+cur_offset(nyyy),' is total res# ',
     &              chainlocal_2_totres(nxxx,nz+cur_offset(nyyy)),
     &              ' in chain# ',nyyy
        enddo
      enddo

C now run scwrl on the structure

      ifirst=1

4444  continue ! come here if scwrl failed the first time

C v3.7 now include option of using "-f scwrl_frame.pdb"
C v4.5 write to scwrl_out.tmp because we'll reorder atoms

C 2022 v1.0 allow option to skip scwrl entirely - for consistency with
C the code we'll try to retain as much as possible of the original
C scwrl-using code, so we'll just copy scwrl_inp.pdb -> scrwl_out.pdb
C here then jump ahead and see if that gets the job done...

      if(skip_scwrl) then
        call system('cp scwrl_inp.pdb scwrl_out.pdb')
        goto 394
      endif

      if(ifirst.eq.1) then
        if(nscwrl_hets.eq.0) then
          sys_string3(1:len_scwrl)=scwrl_path(1:len_scwrl)
          l0=len_scwrl+1
          sys_string3(l0:l0+22)=' -h -t -i scwrl_inp.pdb'
          sys_string3(l0+23:l0+54)=' -s scwrl.fasta -o scwrl_out.tmp'
          sys_string3(l0+55:l0+66)=' > scwrl.out'
        else
          sys_string3(1:len_scwrl)=scwrl_path(1:len_scwrl)
          l0=len_scwrl+1
          sys_string3(l0:l0+22)=' -h -t -i scwrl_inp.pdb'
          sys_string3(l0+23:l0+41)=' -f scwrl_frame.pdb'
          sys_string3(l0+42:l0+73)=' -s scwrl.fasta -o scwrl_out.tmp'
          sys_string3(l0+74:l0+85)=' > scwrl.out'
        endif
      else
        if(nscwrl_hets.eq.0) then
          sys_string3(1:len_scwrl)=scwrl_path(1:len_scwrl)
          l0=len_scwrl+1
          sys_string3(l0:l0+25)=' -h -t -v -i scwrl_inp.pdb'
          sys_string3(l0+26:l0+57)=' -s scwrl.fasta -o scwrl_out.tmp'
          sys_string3(l0+58:l0+69)=' > scwrl.out'
        else
          sys_string3(1:len_scwrl)=scwrl_path(1:len_scwrl)
          l0=len_scwrl+1
          sys_string3(l0:l0+25)=' -h -t -v -i scwrl_inp.pdb'
          sys_string3(l0+26:l0+44)=' -f scwrl_frame.pdb'
          sys_string3(l0+45:l0+76)=' -s scwrl.fasta -o scwrl_out.tmp'
          sys_string3(l0+77:l0+88)=' > scwrl.out'
        endif
      endif
      call system(sys_string3(1:len(trim(sys_string3))))

C get files needed for loopy
C 2022 v1.0 only do this if skip_all_loops is false
C note that since loopy_path points all the way to the executable
C "loopy" we need to trim off the last five characters here and add on
C "model.dir ." instead...

CCC   call system('cp /home/LAB/loopy/test/model.dir .')
      if(.not.skip_all_loops) then
        loopy_copy(1:3)='cp '
        loopy_copy(4:3+len_loopy-5)=loopy_path(1:len_loopy-5)
        loopy_copy(len_loopy-1:len_loopy+9)='model.dir .'
        call system(loopy_copy(1:len_loopy+9))
      endif

C v4.5 now reorder atoms in scwrl_out.pdb so that they are correct!

      greek(1:7)='ABGDEZH'
      char9='XXX X XXX'
      nres_tmp=0
      natm_tmp=0
      open(unit=11,file='scwrl_out.tmp',status='unknown')
      open(unit=12,file='scwrl_out.pdb',status='unknown')
302   read(11,21,end=303)strong

      if(strong(1:4).ne.'ATOM'.and.strong(1:3).ne.'TER') then
        write(12,21)strong
      endif

C note that two conditions trigger when we have reached the end of a
C residue: (1) we read a new res#, or (2) we read a 'TER' statement

      if(strong(18:26).ne.char9.or.strong(1:3).eq.'TER') then
        nres_tmp=nres_tmp+1

C if we just found a new residue then reorder atoms in residue we just read past

        if(nres_tmp.gt.1) then
          idone(1:natm_tmp)=0
          do n2=1,7
            do n1=1,natm_tmp

C we look for the alpha-, beta- etc in the atom name to match - note
C that we also match " " so that we do backbone atoms too - once we've
C matched for an atom in the scwrl_out.tmp file we never write it out
C any more times...

              if(idone(n1).eq.1) cycle
              if(strrng(n1)(15:15).eq.' '.or.
     &           strrng(n1)(15:15).eq.greek(n2:n2)) then
                idone(n1)=1
                write(12,21)strrng(n1)
              endif
            enddo
          enddo
          do n1=1,natm_tmp
            if(idone(n1).eq.0) then
              write(*,*)'no match on scwrl-reorder for: ',strrng(n1)
              stop
            endif
          enddo
        endif

C now, either write out the string originally read if it was 'TER'...

        if(strong(1:3).eq.'TER') then
          write(12,21)strong

C ...or start counting a new residue

        else
          natm_tmp=1
          char9=strong(18:26)
          strrng(natm_tmp)=strong
        endif
      else
        natm_tmp=natm_tmp+1
        strrng(natm_tmp)=strong
      endif
      goto 302
303   close(11)
      close(12)
      write(*,*)
      write(*,*)'done reordering sidechain atoms in scwrl_out.pdb'
      write(*,*)

C 2022 v1.0
C now look for any SG CYS atoms in scwrl_inp.pdb and use them to
C overwrite those in scwrl_out.pdb....
C UPDATE - leave this code in but skip it with a "goto" for now...

      goto 394

      ncys=0
      open(unit=11,file='scwrl_inp.pdb',status='unknown')
1713  read(11,'(a80)',end=1714)char80
      if(char80(1:4).eq.'ATOM') then
        if(char80(14:20).eq.'SG  CYS') ncys=ncys+1
      endif
      goto 1713
1714  rewind(11)

      allocate(xcys(1:ncys))
      allocate(ycys(1:ncys))
      allocate(zcys(1:ncys))

      ncys=0
2713  read(11,'(a80)',end=2714)char80
      if(char80(1:4).eq.'ATOM') then
        if(char80(14:20).eq.'SG  CYS') then
          read(char80,'(30x,3f8.3)')xtmp,ytmp,ztmp
          ncys=ncys+1
          xcys(ncys)=xtmp
          ycys(ncys)=ytmp
          zcys(ncys)=ztmp
        endif
      endif
      goto 2713
2714  close(11)

C now read scwrl_out and replace any CYS SG coords with those from
C scwrl_inp ; after all that, copy the file back to scwrl_out

      ncys=0
      open(unit=11,file='scwrl_out.pdb',status='unknown')
      open(unit=12,file='scwrl_tmp.pdb',status='unknown')
3713  read(11,'(a80)',end=3714)char80
      if(char80(1:4).eq.'ATOM') then
        if(char80(14:20).eq.'SG  CYS') then
          ncys=ncys+1
          write(char24,'(3f8.3)')xcys(ncys),ycys(ncys),zcys(ncys)
          char80(31:54)=char24
        endif
      endif
      write(12,'(a80)')char80
      goto 3713
3714  close(11)
      close(12)

      call system('cp scwrl_tmp.pdb scwrl_out.pdb')

C 2022 v1.0 come here if we used skip_scwrl
C also come here when we skipped doing the CYS coordinate overwrite...

394   continue

C now rename the residues in that output from scwrl...
C note that we also renumber atoms to make sure we don't have
C issues with loopy (scwrl gives -1 for ATOM numbers > 100000)

      natoms_scwrl=0
      nres_scwrl=0
      char9='XXX X XXX'
      open(unit=11,file='scwrl_out.pdb',status='unknown')
      open(unit=12,file='loopy_tmp.pdb',status='unknown')

C v3.2 do a quick read of scwrl_out.pdb to get an estimate of #atoms
C and use this to allocate memory for xs,ys,zs

      nat_tmp=0
3101  read(11,'(a4)',end=3102)char4
      if(char4(1:4).eq.'ATOM'.or.char4(1:4).eq.'HETA') nat_tmp=nat_tmp+1
      goto 3101
3102  rewind(11)

C add num_frags to account for the OXT atoms to be added to each chain
C 2022 v1.0 do this only if skip_scwrl is false
C 2023 v1.0 remove this restriction: do it regardless...

CCC   if(.not.skip_scwrl) then
        nat_tmp=nat_tmp+num_frags
CCC   endif

C 2020 bugfix - make sure allocation only done when ifirst=1...

      if(ifirst.eq.1) then
        allocate(xs(1:nat_tmp))
        allocate(ys(1:nat_tmp))
        allocate(zs(1:nat_tmp))
      endif

C pre-v3.2 code now follows

91    read(11,21,end=92)strong

C TEMP - the TER statements are a complete pain in the ass - they're
C added by scwrl so we'll try to see if we can comment them out here and
C see if that screws things up - it looks like i must want to keep them
C for some reason but i'm not sure what that would be

      if(strong(1:3).eq.'TER') then
CCCCC   write(12,21)strong
        goto 91
      endif

C we also need to skip any SSBOND statements that scwrl adds...

      if(strong(1:6).eq.'SSBOND') goto 91

C otherwise go on and read the atom

      read(strong,93)i,xtmp,ytmp,ztmp
93    format(22x,i4,4x,3f8.3)

C v2.0 - here we try to make things work when residue numbers exceed
C 10,000 - previously we used to just read the residue number direct
C from scwrl_out.pdb and then convert that into an actual total residue
C number - that doesn't seem to work - so now we'll just keep track of
C every new residue we find in scwrl_out.pdb and use that as the counter
C - this should be easy to do as scwrl doesn't add "new" residues

      if(strong(18:26).ne.char9) then
        nres_scwrl=nres_scwrl+1
        char9=strong(18:26)
      endif

ccc   ijk=my_tot_res(i)
      ijk=my_tot_res(nres_scwrl)

C we need to make sure that we skip residues that were in the
C scwrl_out.pdb but that are going to be part of loops - we don't want
C them to already be present when we try to build loops or tails...
C note that it's somewhat incongruous to have these residues included at
C the scwrl stage but then deleted later, but maybe it doesn't matter...

C note that i_done_res used to be cosmetic, but now it's used to
C differentiate residues that we want to rebuild but for which we
C already have coordinates...

      if(i_am_original(ijk).eq.0) then
        if(i_done_res(ijk).eq.0) then
          if(strong(13:16).eq.' N  ') ! write 1st atom to screen
     &    write(*,8877)i,ijk
8877      format('skipping transfer of scwrl-res ',i6,
     &           ' i.e. real-res ',i6,' to loopy_tmp.pdb')
          goto 91
        endif
      endif

C v2.1 - I don't believe that we need to change the res# when going from
C scwrl to loopy now so comment the following lines out:

C     write(num_of_res,'(i4)')mod(ijk,10000)
C     strong(23:26)=num_of_res(1:4)

C v2.1 - the following must be a defunct line from *ages* ago...
ccc   strong(22:22)=' ' ! need to remove chain ID for loopy...

C check to see if we're on an actual C-terminus and if so keep coords of
C the atoms that we'll need to place the OXT atom... note that this is
C necessary because scwrl has been set to ignore termini (which is the
C correct decision on balance) and this means it ignores the real OXT
C atoms when it finds them (and we never add on OXT atoms for those
C C-terminal residues that are non-terminal *in the template*). The
C following lines, and the write statements a page down, fix this issue.

C 2022 v1.0 only do this if skip_scwrl is false
C 2023 v1.0 do this regardless for now - otherwise we lose the CTER OXT

C     if(.not.skip_scwrl) then
        if(strong(13:16).eq.' O  ') then
          do i=1,num_act_cters
            if(ijk.eq.ipos_act_cters(i)) then
              a(1)=xtmp
              a(2)=ytmp
              a(3)=ztmp
            endif
          enddo
        elseif(strong(13:16).eq.' CA ') then
          do i=1,num_act_cters
            if(ijk.eq.ipos_act_cters(i)) then
              b(1)=xtmp
              b(2)=ytmp
              b(3)=ztmp
            endif
          enddo
        elseif(strong(13:16).eq.' C  ') then
          do i=1,num_act_cters
            if(ijk.eq.ipos_act_cters(i)) then
              c(1)=xtmp
              c(2)=ytmp
              c(3)=ztmp
            endif
          enddo
        endif
C     endif

      natoms_scwrl=natoms_scwrl+1
c     write(12,21)strong
      write(12,'("ATOM ",i6,a69)')natoms_scwrl,strong(12:80)

C store the coordinates as well for subsequent trimming of HETATMs

      xs(natoms_scwrl)=xtmp
      ys(natoms_scwrl)=ytmp
      zs(natoms_scwrl)=ztmp

C if we just wrote out an O atom we can check if it was a C-terminus and
C if necessary write out the OXT atom to go with it also...
C note that we will already know a(1:3),b(1:3),c(1:3) by this point too

C 2022 v1.0 only do this if skip_scwrl is false
C 2023 v1.0 do this regardless for now - otherwise we lose the CTER OXT

C     if(.not.skip_scwrl) then
        if(strong(13:16).eq.' O  ') then
          do i=1,num_act_cters
            if(ijk.eq.ipos_act_cters(i)) then
              bondcur=1.25
              anglcur=118.0*3.14159265/180.0
              dihecur=180.0*3.14159265/180.0
              st=sin(anglcur)
              zpd=-bondcur*cos(anglcur)
              xpd= bondcur*cos(dihecur)*st
              ypd= bondcur*sin(dihecur)*st
              do l=1,3
                vca(l)=a(l)-c(l)
                vcb(l)=b(l)-c(l)
              enddo
              call acrosb(vca,vcb,yp)
              call acrosb(vcb,yp,xp)
              call acrosb(xp,yp,zp)
              xoxt=xp(1)*xpd+yp(1)*ypd+zp(1)*zpd+c(1)
              yoxt=xp(2)*xpd+yp(2)*ypd+zp(2)*zpd+c(2)
              zoxt=xp(3)*xpd+yp(3)*ypd+zp(3)*zpd+c(3)
              strong(13:16)=' OXT'
              natoms_scwrl=natoms_scwrl+1
              write(12,'("ATOM ",i6,a19,3f8.3,a26)')natoms_scwrl,
     &          strong(12:30),xoxt,yoxt,zoxt,strong(55:80)

C v4.1 BUGFIX - we should store xs,ys,zs for this new atom also
C this was not stored previously but probably wasn't being used anyway
C as we had skipped the xs-based attempt to identify HETATM entries to
C keep - now that we have reinstated the latter we need to make sure
C that the xs,ys,zs coords are correctly stored

              xs(natoms_scwrl)=xoxt
              ys(natoms_scwrl)=yoxt
              zs(natoms_scwrl)=zoxt

            endif
          enddo
        endif
C     endif

      goto 91
92    close(11)

C at this stage we need to put in a check to see if scwrl worked...

      if(natoms_scwrl.eq.0.and.ifirst.eq.1) then
        write(*,*)
        write(*,*)'WARNING WARNING WARNING'
        write(*,*)'scwrl seems to have generated no atoms'
        write(*,*)'we will try to run again using the "-v" option'
        write(*,*)'note that this may be a sign of a serious problem!'
        write(*,*)
        ifirst=2
        goto 4444
      elseif(natoms_scwrl.eq.0.and.ifirst.eq.2) then
        write(*,*)
        write(*,*)'scwrl seems to have generated no atoms'
        write(*,*)'even when using the "-v" option'
        write(*,*)'you should look in scwrl.out to debug this'
        write(*,*)'quitting now :('
        write(*,*)
        stop
      endif

C v2.0 we also want to check that we counted residues correctly...

      nres_scwrl_out=nres_scwrl
      if(nres_scwrl_inp.ne.nres_scwrl_out) then
        write(*,*)
        write(*,*)'nres_scwrl_inp = ',nres_scwrl_inp
        write(*,*)'nres_scwrl_out = ',nres_scwrl_out
        write(*,*)'these numbers should be identical - you may be '
        write(*,*)'running into problems due to having >10,000 res'
        write(*,*)'quitting now :('
        write(*,*)
        stop
      endif

C here we can add a check to trim the list of HETATM entries to include
C only those that are within 4A of an atom in scwrl_out.pdb - or that
C are within 4A of an atom that is...

C v3.1 skip this step for now - TEMP TEMP TEMP
C v4.1 this really should be included so we'll speed it up now by
C putting scwrl atoms on to a grid - note that we don't use a grid for
C the already-selected hetatms - we just do a regular loop
C UPDATE - let's make this optional using iskiphetcheck...

      if(iskiphetcheck.eq.1) then
        write(*,*)
        write(*,*)'WARNING - skipping the HETATM distance check'
        write(*,*)'we will include *ALL* HETATMs for now...'
        write(*,*)

C v4.5 9 April - we can skip the iterative scheme and still select all
C HETATMs correctly if we set ih() and ih_original() both to 1 here

        do n=1,num_het_atms
          ih(n)=1
          ih_original(n)=1
        enddo

        goto 4994
      endif

      xmin= 999999.9
      ymin= 999999.9
      zmin= 999999.9
      xmax=-999999.9
      ymax=-999999.9
      zmax=-999999.9
      do m=1,natoms_scwrl
        xmin=min(xmin,xs(m))
        ymin=min(ymin,ys(m))
        zmin=min(zmin,zs(m))
        xmax=max(xmax,xs(m))
        ymax=max(ymax,ys(m))
        zmax=max(zmax,zs(m))
      enddo
      xlen=xmax-xmin+10.0
      ylen=ymax-ymin+10.0
      zlen=zmax-zmin+10.0
      xinv=1.0/4.0 
      yinv=1.0/4.0 
      zinv=1.0/4.0 
      nx=int(xlen/4.0)+1
      ny=int(ylen/4.0)+1
      nz=int(zlen/4.0)+1
      allocate(temp_grid(1:nx,1:ny,1:nz))
      temp_grid=0
      ninbox_max=0

      write(*,*)'xmin etc ',xmin,ymin,zmin
      write(*,*)'xmax etc ',xmax,ymax,zmax
      write(*,*)'xlen etc ',xlen,ylen,zlen
      write(*,*)'xinv etc ',xinv,yinv,zinv
      write(*,*)'nx   etc ',nx,ny,nz

      do m=1,natoms_scwrl
        xt=xs(m)
        yt=ys(m)
        zt=zs(m)
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        write(*,*)'check ',m,i1,j1,k1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              ninbox_max=max(ninbox_max,temp_grid(i1,j1,k1))
            enddo
          enddo
        enddo
      enddo
 
      write(*,*)'for HETATM check, ninbox_max = ',ninbox_max
      write(*,*)'allocating memory for grid'

      allocate(xatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))
      allocate(yatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))
      allocate(zatm_grid(1:ninbox_max,1:nx,1:ny,1:nz))

      temp_grid=0

      do m=1,natoms_scwrl
        xt=xs(m)
        yt=ys(m)
        zt=zs(m)
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              xatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=xt
              yatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=yt
              zatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=zt
            enddo
          enddo
        enddo
      enddo

      write(*,*)'done putting everyone on the grid'

      num_het_atms_final=0
      ih=0
      do n=1,num_het_atms

        xt=xh(n)
        yt=yh(n)
        zt=zh(n)
        iwant=0
        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1

C HETATMs might be off the grid (which only fits scwrl_out.pdb atoms)
C so be prepared to skip this HETATM...

        if(i1.lt.1) cycle
        if(j1.lt.1) cycle
        if(k1.lt.1) cycle
        if(i1.gt.nx) cycle
        if(j1.gt.ny) cycle
        if(k1.gt.nz) cycle

        do l1=1,temp_grid(i1,j1,k1)
          xu=xatm_grid(l1,i1,j1,k1) 
          yu=yatm_grid(l1,i1,j1,k1) 
          zu=zatm_grid(l1,i1,j1,k1) 
          dist=(xt-xu)**2+
     &         (yt-yu)**2+
     &         (zt-zu)**2
          if(dist.le.16.0) then
            iwant=1
            exit
          endif
        enddo

        if(iwant.eq.1) then
          ih(n)=1 ! make sure we don't do this HETATM again
          num_het_atms_final=num_het_atms_final+1
ccc       si(num_het_atms_final)=sh(n) ! v4.5
ccc       mi(num_het_atms_final)=mh(n) ! v4.5
ccc       xi(num_het_atms_final)=xh(n)
ccc       yi(num_het_atms_final)=yh(n)
ccc       zi(num_het_atms_final)=zh(n)
        endif

      enddo

      deallocate(temp_grid)
      deallocate(xatm_grid)
      deallocate(yatm_grid)
      deallocate(zatm_grid)

      write(*,*)'after looking at protein atoms only '
      write(*,*)'current # of HETATMs to include = ',
     &           num_het_atms_final

C now we need to keep iterating and look for additional HETATMs that
C should be added to the list due to proximity to others in the list

5052  continue

      num_het_atms_previous=num_het_atms_final

      do n=1,num_het_atms

        if(ih(n).eq.1) cycle ! skip if we already did this one
        iwant=0

ccc     do m=1,num_het_atms_final
        do m=1,num_het_atms

          if(ih(m).eq.0) cycle ! skip if we didn't already find it!
                               ! we only want to compare atom n with
                               ! atoms that we've already selected

c         dist=(xh(n)-xi(m))**2+
c    &         (yh(n)-yi(m))**2+
c    &         (zh(n)-zi(m))**2
          dist=(xh(n)-xh(m))**2+
     &         (yh(n)-yh(m))**2+
     &         (zh(n)-zh(m))**2
          if(dist.le.16.0) then
            iwant=1
            goto 5053
          endif

        enddo
5053    continue

        if(iwant.eq.1) then
          ih(n)=1 ! make sure we don't do this HETATM again
          num_het_atms_final=num_het_atms_final+1
ccc       si(num_het_atms_final)=sh(n) ! v4.5
ccc       mi(num_het_atms_final)=mh(n) ! v4.5
ccc       mi(num_het_atms_final)=n     ! v4.5 store original het#
ccc       xi(num_het_atms_final)=xh(n)
ccc       yi(num_het_atms_final)=yh(n)
ccc       zi(num_het_atms_final)=zh(n)
        endif

      enddo
      
      if(num_het_atms_final.gt.num_het_atms_previous) then
        write(*,*)'current # of HETATMs to include = ',
     &             num_het_atms_final
        goto 5052
      else
        write(*,*)'FINAL   # of HETATMs to include = ',
     &             num_het_atms_final

C store ih_original for all HETATM entries - we'll use this to reset ih
C after we've built each MODEL

        do n=1,num_het_atms
          ih_original(n)=ih(n)
        enddo

      endif

C if we got here then we've finished making the final list of HETATMs -
C we can now copy these ones back to the final list...

      ntmp=0
      open(unit=21,file='HETERO_original.pdb',status='unknown')
      open(unit=29,file='HETERO_temporary.pdb',status='unknown')
5054  read(21,'(a80)',end=5055)strong
      ntmp=ntmp+1
      if(ih(ntmp).eq.1) write(29,'(a80)')strong
      goto 5054
5055  close(21)
      close(29)
      
      call system('mv HETERO_temporary.pdb HETERO_original.pdb')

C if we didn't ever set ih=1 for a HETATM then we need to set mh=0 for
C that HETATM so that we don't write it out at the very end...
C UPDATE just use ih() instead

ccc   do n=1,num_het_atms
ccc     if(ih(n).eq.0) mh(n)=0
ccc   enddo
ccc
ccc   num_het_atms=0
ccc   do n=1,num_het_atms_final
ccc     num_het_atms=num_het_atms+1
ccc     sh(n)=si(n) ! v4.5
ccc     mh(n)=mi(n) ! this will be the original hetatm # of this one
ccc     xh(n)=xi(n)
ccc     yh(n)=yi(n)
ccc     zh(n)=zi(n)
ccc   enddo

4994  continue

C before closing loopy_tmp.pdb we need to add on the HET atom entries 
C that weren't deleted as GLY residues - this is an easy fix to make 
C sure that loops are built in the knowledge of other groups present

C v4.1 here is another example of an old complaint that is now fixed:

CCC note that I *wish* we could do the same with the scwrl calculation but
CCC there doesn't seem to be an easy way to do it :(
CCC UPDATE - I think there is - it's just not obvious

C v3.1 note that the chainID and res#s assigned to the fake GLY residues
C should not matter ultimately as they'll be filtered out later

C v4.7 the find_stretches code fails if the HETERO file has asterisks
C for atom numbers - so here we'll explicitly write out atom numbers
C that cannot possibly exceed 100000 - formerly we were not adjusting
C columns 7-12 at all...

      ntmp=0
      nhetres=-999
      ngly=0
      open(unit=21,file='HETERO_original.pdb',status='unknown')
71    read(21,21,end=73)strong
      ntmp=ntmp+1
      if(ntmp.eq.5) then ! deal with atoms in chunks of 4
        nhetres=nhetres+1
        ntmp=1
      endif

      strong(1:6)='ATOM  '
      write(num_of_res,'(i4)')mod(nhetres,10000)
      ijk=mod(num_frags+1,chainIDmax)
      if(ijk.eq.0) ijk=chainIDmax
      strong(22:22)=blah(ijk:ijk) ! this adds one to chain name
      strong(23:26)=num_of_res
      strong(55:60)='-99.99'      ! use this later to delete these

      if(ntmp.eq.1) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' N   GLY',strong(21:80)
      elseif(ntmp.eq.2) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' CA  GLY',strong(21:80)
      elseif(ntmp.eq.3) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' C   GLY',strong(21:80)
      elseif(ntmp.eq.4) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' O   GLY',strong(21:80)
      endif
72    format(a12,a8,a60)
      goto 71
73    close(21)

C we need to make sure we have a whole GLY residue at the end...

      if(ntmp.eq.1) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' CA  GLY',strong(21:80)
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' C   GLY',strong(21:80)
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' O   GLY',strong(21:80)
      elseif(ntmp.eq.2) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' C   GLY',strong(21:80)
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' O   GLY',strong(21:80)
      elseif(ntmp.eq.3) then
        ngly=ngly+1
        nwrt=mod(ngly,100000)
        write(strong(7:11),'(i5)')nwrt
        write(12,72)strong(1:12),' O   GLY',strong(21:80)
      endif
 
C at this point we've added on the fake GLY residues for HET atoms so we
C can close loopy_tmp.pdb and proceed...

      write(*,*)'# fake GLY atoms added to loopy_tmp.pdb = ',ngly
 
      close(12)     

C v2.1 keep a copy of loopy_tmp.pdb

ccc   call system('cp loopy_tmp.pdb loopy_tmp.pdb.original')
      open(unit=11,file='loopy_tmp.pdb',status='unknown') ! input
      open(unit=12,file='loopy_tmp.pdb.original',status='unknown') ! output
8581  read(11,'(a80)',end=8582)char80
      write(12,'(a80)')char80
      goto 8581
8582  close(11)
      close(12)

C v3.6 - option here to skip building all loops - make sure that
C loopy_tmp.pdb is copied to loopy_inp.pdb as this will be the final
C output file generated during loop building...

      if(skip_all_loops) then

ccc     call system('cp loopy_tmp.pdb loopy_inp.pdb')
        open(unit=11,file='loopy_tmp.pdb',status='unknown') ! input
        open(unit=12,file='loopy_inp.pdb',status='unknown') ! output
8571    read(11,'(a80)',end=8572)char80
        write(12,'(a80)')char80
        goto 8571
8572    close(11)
        close(12)

        goto 7654

      endif

C now build in all the missing loops one at a time...

      idoneallocating=0
      do i=1,num_tot_loops

ccc     if(i.eq.21) stop ! TEMP TEMP TEMP

        write(*,*)

C immediately skip if we already know that this is dead loop

        if(ibeg_loc_loops(i).eq.-999) then
          write(*,*)'skipping merged loop # ',i
          cycle ! skip merged loops
        endif

C v2.1 keep the fragment# of this loop in a temporary variable

        my_frag=ifrg_tot_loops(i)

C v1.9 here is where we can read loopy_tmp.pdb and assess if a loop is
C even possible from a geometric point of view: if the distance between
C the beg and end of the loop exceeds plausibility then we could decide
C to just skip it - useful for handling multi-template cases...

        loopstring(1:14)='rm strch.txt; '
        loopstring(15:15+len_code-1)=code_path(1:len_code)
        l0=15+len_code
        loopstring(l0:l0+35)='/find_stretches_loopy_pdb_file.exe  '
        loopstring(l0+36:l0+68)='loopy_tmp.pdb          strch.txt '
        loopstring(l0+69:l0+69)=ichn_tot_loops(i)
C v2.1  write(a1,'(i5)')ibeg_tot_loops(i)-1
        write(a1,'(i5)')ibeg_loc_loops(i)+cur_offset(my_frag)-1
        j1=len(trim(adjustl(a1)))
C v2.1  write(a2,'(i5)')iend_tot_loops(i)+1
        write(a2,'(i5)')iend_loc_loops(i)+cur_offset(my_frag)+1
        j2=len(trim(adjustl(a2)))
        loopstring(l0+70:l0+70)=' '
C v2.1  write(a1,'(i5)')ibeg_tot_loops(i)-1
        write(a1,'(i5)')ibeg_loc_loops(i)+cur_offset(my_frag)-1
        loopstring(l0+71:l0+75)=a1
        loopstring(l0+76:l0+76)=' '
C v2.1  write(a1,'(i5)')iend_tot_loops(i)+1
        write(a1,'(i5)')iend_loc_loops(i)+cur_offset(my_frag)+1
        loopstring(l0+77:l0+81)=a1
        loopstring(l0+82:l0+82)=' '
        write(a1,'("     ")')
        loopstring(l0+83:l0+87)=a1
        write(*,*)"system call: ",loopstring(1:l0+87)
        call system(loopstring(1:l0+87))

C IS THE ABOVE OR BELOW CODE ERROR-PRONE? I'M SUSPICIOUS THAT SOMETIMES
C STRCH.TXT MAY NOT HAVE STUFF IN IT WHEN IT SHOULD...
C UPDATE - THIS MAY BE BECAUSE OF A DIFFERENT BUG RELATING TO RESIDUES
C BEING WRITTEN OUT MULTIPLE TIMES IF PRESENT IN ORIGINAL TEMPLATE AND
C BEING BUILT BY LOOPY SIMULTANEOUSLY...

C UPDATE - it may have been error-prone with very large systems due to
C array allocation issues - the code is now allocated to 50 million 

        call system('touch strch.txt ; wc -l strch.txt > strch.tmp')
        open(unit=81,file='strch.tmp',status='unknown')
        read(81,*)ijklm
        if(ijklm.ne.0) then
          ibeg_loc_loops(i)=-999
          write(*,*)'loop has to stretch insane distance - ignore it'
        endif
        close(81)

606     continue ! come back here to have another go at a loop

C check again to make sure this hasn't become a dead loop recently

        if(ibeg_loc_loops(i).eq.-999) then
          write(*,*)'skipping merged loop # ',i
          cycle ! skip merged loops
        endif

C v2.1 need to make loopy_inp.pdb from loopy_tmp.pdb - just include
C those residues that are "near" the possible loop...
C v2.5 but exclude any that are part of the loop that we're building - I
C don't seem to have found it necessary to do this in earlier versions
C but we have a clear case where it goes crazy...
C note that when the new version of loopy_tmp.pdb is made at the end of
C this iteration it will naturally omit the extra copies...

C v3.1 to solve problems that we seem to have with loopy when residue
C numbers are in the thousands we will change the chainID and the
C residue numbers of all atoms that get put in loopy_inp.pdb - the only
C ones that we really need to keep track of are the ones of the loop
C that we want to make, so we'll record those...

        open(unit=11,file='loopy_tmp.pdb',status='unknown')
        open(unit=12,file='loopy_inp.pdb',status='unknown')

C first, we read the CA coords of the residues -1 and +1 of the loop
C that we're building - these will be used to define an axis-aligned
C bounding box (AABB) that we'll use to decide which residues to keep

        jbeg=ibeg_loc_loops(i)+cur_offset(my_frag)-1
        jend=iend_loc_loops(i)+cur_offset(my_frag)+1
9001    read(11,'(a66)',end=9002)string66
        read(string66,'(21x,a,i4)')chainID,jkl
        if(chainID.eq.ichn_tot_loops(i).and.jkl.eq.jbeg.and.
     &     string66(13:16).eq.' CA ') then
          read(string66,'(30x,3f8.3)')xq(1),yq(1),zq(1)
        elseif(chainID.eq.ichn_tot_loops(i).and.jkl.eq.jend.and.
     &     string66(13:16).eq.' CA ') then
          read(string66,'(30x,3f8.3)')xq(2),yq(2),zq(2)
          goto 9002 ! we know we're done reading if we found jend
        endif
        goto 9001
9002    rewind(11)
        xmi=min(xq(1),xq(2))
        ymi=min(yq(1),yq(2))
        zmi=min(zq(1),zq(2))
        xma=max(xq(1),xq(2))
        yma=max(yq(1),yq(2))
        zma=max(zq(1),zq(2))
        xmi=xmi-20.0
        ymi=ymi-20.0
        zmi=zmi-20.0
        xma=xma+20.0
        yma=yma+20.0
        zma=zma+20.0

        write(*,*)
        write(*,*)'for loop# ',i,' out of ',num_tot_loops
        write(*,*)'in fragment# ',my_frag
        write(*,*)
        write(*,*)'keeping atoms from loopy_tmp.pdb within AABB '
        write(*,*)'xmi & xma ',xmi,xma
        write(*,*)'ymi & yma ',ymi,yma
        write(*,*)'zmi & zma ',zmi,zma
        write(*,*)

C now that we know where the axis-aligned bounding box is for the loop
C we can go ahead and find all residues in loopy_tmp.pdb that have at
C least one atom within the AABB...
C note that we keep a list here of all unique combinations of chainID
C and residue number - they are held as "cwithin" and "rwithin", resp.
C we then go and find all atoms that share these attributes - this
C allows us to ensure that all residues given to loopy are complete...

C v2.5 exclude any residues that are part of the loop to be built...
C for this purpose, redefine jbeg and jend here...

        jbeg=ibeg_loc_loops(i)+cur_offset(my_frag)
        jend=iend_loc_loops(i)+cur_offset(my_frag)

C v3.1 for writing to loopy_inp.pdb we will make the chainID of the loop
C be "A" and the chainID of all other residues be "B" - since we're
C using an AABB we shouldn't need more than 999 residues in B...
C also change jbeg to jbeg_loopy and jend to jend_loopy

        jbeg_loopy=ibeg_loc_loops(i)
        jend_loopy=iend_loc_loops(i)

        nwithin=0
        natoms_tmp=0
        natoms_gly=0 ! #fake GLY atoms that we want
        nres_gly=0   ! #fake GLY residues that we want

        if(idoneallocating.eq.0) then
          allocate(iwant_this_gly(1:ngly))
          idoneallocating=1
        endif

        iwant_this_gly=0 ! set this to zero for all fake GLY atoms

9003    read(11,'(a66)',end=9004)string66
        read(string66,'(21x,a,i4,4x,3f8.3,a6)')chainID,jkl,x,y,z,char6

C increment both counters before possibly skipping this line

        natoms_tmp=natoms_tmp+1
        if(char6.eq.'-99.99') natoms_gly=natoms_gly+1

C v2.5 here's where we skip those intra-loop residues

        if(chainID.eq.ichn_tot_loops(i).and.
     &     jkl.ge.jbeg.and.
     &     jkl.le.jend) goto 9003

        if(x.ge.xmi.and.x.le.xma.and.y.ge.ymi.and.y.le.yma.and.
     &     z.ge.zmi.and.z.le.zma) then

C check the occupancy - if it's -99.99 then this is a fake GLY residue
C if so, then we want to keep this atom and its 3 partner atoms in the
C same fake GLY residue - note that these are always written in groups
C of 4 (see above) so there is no question about their numbering

          if(char6.eq.'-99.99') then 
            read(string66,'(12x,a4)')char4
            if(char4.eq.' N  ') then

C if it's N then it's the first atom of this fake residue so we need to
C increment nres_gly and assign it to all 4 atoms of this residues

              nres_gly=nres_gly+1
              iwant_this_gly(natoms_gly)  =nres_gly
              iwant_this_gly(natoms_gly+1)=nres_gly
              iwant_this_gly(natoms_gly+2)=nres_gly
              iwant_this_gly(natoms_gly+3)=nres_gly
            elseif(char4.eq.' CA ') then

C if it's CA then we should first check if the first atom of this
C residue has been assigned a non-zero value - if it has then we're
C done; otherwise, we need to increment nres_gly and assign it to all 4
C atoms of this fake residue

              if(iwant_this_gly(natoms_gly-1).ne.0) goto 9003
              nres_gly=nres_gly+1
              iwant_this_gly(natoms_gly-1)=nres_gly
              iwant_this_gly(natoms_gly)  =nres_gly
              iwant_this_gly(natoms_gly+1)=nres_gly
              iwant_this_gly(natoms_gly+2)=nres_gly
            elseif(char4.eq.' C  ') then

C logic here is the same as above

              if(iwant_this_gly(natoms_gly-2).ne.0) goto 9003
              nres_gly=nres_gly+1
              iwant_this_gly(natoms_gly-2)=nres_gly
              iwant_this_gly(natoms_gly-1)=nres_gly
              iwant_this_gly(natoms_gly)  =nres_gly
              iwant_this_gly(natoms_gly+1)=nres_gly
            elseif(char4.eq.' CA ') then

C logic here is the same as above

              if(iwant_this_gly(natoms_gly-3).ne.0) goto 9003
              nres_gly=nres_gly+1
              iwant_this_gly(natoms_gly-3)=nres_gly
              iwant_this_gly(natoms_gly-2)=nres_gly
              iwant_this_gly(natoms_gly-1)=nres_gly
              iwant_this_gly(natoms_gly)  =nres_gly
            endif 
            goto 9003
          endif

C if we get here then it's a regular atom - here we decide whether it's
C a new chain/res combination or not

          iamnew=1
          do n8=1,nwithin
            if(chainID.eq.cwithin(n8).and.jkl.eq.rwithin(n8)) then
              iamnew=0
              exit
            endif
          enddo
          if(iamnew.eq.1) then
            nwithin=nwithin+1
            cwithin(nwithin)=chainID
            rwithin(nwithin)=jkl
            lwithin(nwithin)=nwithin ! new res# for loopy_inp.pdb
          endif
        endif

        goto 9003
9004    rewind(11)

        write(*,*)'#res to include in loopy_inp.pdb= ',nwithin
        write(*,*)'#fake GLY res to include in loopy_inp.pdb= ',nres_gly

        natoms_tmp=0
        natoms_gly=0 ! #fake GLY atoms that we want
        natoms_out=0 ! loopy needs small atom #s so renumber

        n2write=0
9005    read(11,'(a66)',end=9006)string66
        read(string66,'(21x,a,i4,4x,3f8.3,a6)')chainID,jkl,x,y,z,char6

C increment both counters before possibly skipping this line

        natoms_tmp=natoms_tmp+1
        if(char6.eq.'-99.99') natoms_gly=natoms_gly+1

C if this atom is in a fake GLY residue that we want then it is given
C the appropriate res# and given a chainID of "D"

        if(char6.eq.'-99.99') then
          if(iwant_this_gly(natoms_gly).ne.0) then
            natoms_out=natoms_out+1
            write(char6,'(i6)')natoms_out
            string66(6:11)=char6
            string66(22:22)="D"
            write(char4,'(i4)')iwant_this_gly(natoms_gly)
            string66(23:26)=char4
            write(12,'(a66)')string66
            n2write=n2write+1
            goto 9005
          endif
        endif

        do n8=1,nwithin
          if(chainID.eq.cwithin(n8).and.jkl.eq.rwithin(n8)) then

C v3.1 change the chainID to A for residues in same chain as loop;
C change the res# for residues in the same chain as loop; change the
C chainID to B for all other residues with the same chainID (but not
C actually in the smae chain); change the res# for all other 
C residues to that determined above...
C the only things we need to be careful about converting later are the
C residue numbers in the same chain as the loop...

            if(chainID.eq.ichn_tot_loops(i).and.
     &         jkl.ge.cur_offset(my_frag)+1.and.
     &         jkl.le.cur_offset(my_frag)+num_res_in_frag(my_frag)) then
              string66(22:22)="A"
              jkl_tmp=jkl-(jbeg-jbeg_loopy)
              write(char4,'(i4)')jkl_tmp
              string66(23:26)=char4
            elseif(chainID.eq.ichn_tot_loops(i).and.
     &       (jkl.le.cur_offset(my_frag).or.
     &        jkl.gt.cur_offset(my_frag)+num_res_in_frag(my_frag))) then
              string66(22:22)="B"
              write(char4,'(i4)')lwithin(n8)
              string66(23:26)=char4
            else
              string66(22:22)="C"
              write(char4,'(i4)')lwithin(n8)
              string66(23:26)=char4
            endif

            natoms_out=natoms_out+1
            write(char6,'(i6)')natoms_out
            string66(6:11)=char6
            write(12,'(a66)')string66
            n2write=n2write+1
            goto 9005
          endif
        enddo
        goto 9005
9006    close(11)
        close(12)

        write(*,*)'#atm to include in loopy_inp.pdb= ',n2write

C v2.5 also keep a copy of every loop immediately befor making...
C v4.5 - not yet implemented replacement of system call

        sys_stringw(1:37)='cp loopy_inp.pdb loopy_SUBMIT_XXX.pdb'
        write(char3,'(i3)')i
        if(char3(1:1).eq.' ') char3(1:1)='0'
        if(char3(2:2).eq.' ') char3(2:2)='0'
        sys_stringw(31:33)=char3
        call system(sys_stringw(1:37))

        do lll=1,1000
          sys_string4(lll:lll)=' '
        enddo
C v3.1  write(a1,'(i5)')ibeg_loc_loops(i)+cur_offset(my_frag)
        write(a1,'(i5)')jbeg_loopy
        j1=len(trim(adjustl(a1)))
C v3.1  write(a2,'(i5)')iend_loc_loops(i)+cur_offset(my_frag)
        write(a2,'(i5)')jend_loopy
        j2=len(trim(adjustl(a2)))
        write(a3,'(i5)')nloopy_confs ! # of loop confs to build
        j3=len(trim(adjustl(a3)))
        sys_string4(1:12)='timeout  1h '
        sys_string4(13:13+len_loopy-1)=loopy_path(1:len_loopy)
        sys_string4(13+len_loopy:13+len_loopy+3)=' -o='
        l0=len_loopy+17
        sys_string4(l0:l0+j1-1)=trim(adjustl(a1))
        sys_string4(l0+j1:l0+j1)='-'
        sys_string4(l0+j1+1:l0+j1+j2)=trim(adjustl(a2))
        sys_string4(l0+j1+j2+1:l0+j1+j2+17)=' -t=b    -s=c -r='
        k1=ibeg_tot_loops(i)
        k2=iend_tot_loops(i)
        j1=l0+j1+j2+18
        j2=j1+iend_tot_loops(i)-ibeg_tot_loops(i)
        sys_string4(j1:j2)=seq_tot_res(k1:k2)
        sys_string4(j2+1:j2+18)=' loopy_inp.pdb -c='
C v3.1  sys_string4(j2+19:j2+19)=ichn_tot_loops(i)
        sys_string4(j2+19:j2+19)="A"
        sys_string4(j2+20:j2+23)=' -n='
        sys_string4(j2+24:j2+24+j3-1)=trim(adjustl(a3))
        sys_string4(j2+24+j3:j2+36+j3)=' > temp1 2>&1'
        k1=len(trim(sys_string4))
        write(*,*)'CALLING ',sys_string4(1:k1)
        call system(sys_string4(1:k1))

C update i_am_by_loopy - note that these assignments will be made over
C and over again if loopy is called repeatedly but this is not a problem
 
        do ijk=ibeg_tot_loops(i),iend_tot_loops(i)
          i_am_by_loopy(ijk)=1
        enddo

C start what is in effect a loop over the conformations made by loopy -
C we'll keep looking through them till we find one that's acceptable

        imad(1:10)=0
        jmad(1:10)=0
        do nmad=1,10
8765      call random_number(u1)
          iloop=int(u1*10)+1
          if(imad(iloop).eq.0) then
            jmad(nmad)=iloop
            imad(iloop)=1
          else
            goto 8765
          endif
          write(*,*)'check randomized loop selection ',nmad,jmad(nmad)
        enddo

        nmade=0
9876    continue   ! come back here to check next conf made by loopy
        nmade=nmade+1 ! stores the conf# that we're checking
                      ! in v1.2 this just tracks how many confs we've
                      ! checked to see if they work...

        nmad2=jmad(nmade) ! convert the randomized loop we want

C read the output pdb file from loopy and make sure that the loop is
C actually a complete one

        iloopydone=1 ! assume that loopy worked well for this conf

C take steps to verify whether loopy_inp_looper_.pdb got made

        sys_string6(1:28)='touch loopy_inp_looper_X.pdb'
        write(sys_string6(24:24),'(i1)')nmad2-1
        call system(sys_string6(1:28))

        sys_string8(1:36)='wc -l loopy_inp_looper_X.pdb > temp3'
        write(sys_string8(24:24),'(i1)')nmad2-1
        call system(sys_string8(1:36))

C we read the result of the "wc -l " command - if it's zero then the
C loopy...pdb file was not made at all...

        open(unit=12,file='temp3',status='unknown')
        read(12,*)iloopydone ! =0 if the file was not made
        close(12)

C assuming that iloopydone<>0 we can check for clashes...

        loopstring(1:14)='rm clash.txt; '
        loopstring(15:15+len_code-1)=code_path(1:len_code)
        l0=15+len_code
        loopstring(l0:l0+35)='/find_clashes_in_loopy_pdb_file.exe '
        loopstring(l0+36:l0+69)='loopy_inp_looper_X.pdb clash.txt  '
        write(loopstring(l0+53:l0+53),'(i1)')nmad2-1
C v3.1  loopstring(111:111)=ichn_tot_loops(i)
        loopstring(l0+70:l0+70)="A"
C v3.1  write(a1,'(i5)')ibeg_loc_loops(i)+cur_offset(my_frag)
        write(a1,'(i5)')jbeg_loopy
        j1=len(trim(adjustl(a1)))
C v3.1  write(a2,'(i5)')iend_loc_loops(i)+cur_offset(my_frag)
        write(a2,'(i5)')jend_loopy
        j2=len(trim(adjustl(a2)))
        loopstring(l0+71:l0+71)=' '
C v3.1  write(a1,'(i5)')ibeg_loc_loops(i)+cur_offset(my_frag)
        write(a1,'(i5)')jbeg_loopy
        loopstring(l0+72:l0+76)=a1
        loopstring(l0+77:l0+77)=' '
C v3.1  write(a1,'(i5)')iend_loc_loops(i)+cur_offset(my_frag)
        write(a1,'(i5)')jend_loopy
        loopstring(l0+78:l0+82)=a1
        loopstring(l0+83:l0+83)=' '
C v4.7 BUGFIX 
CCC     write(a1,'(f5.2)')clash_distance
CCC     write(a1,'(f5.2)')clash_intra
C 2023_v1.1 - call with max of clash_dist1 & clash_dist2
        write(a1,'(f5.2)')clash_max
        loopstring(l0+84:l0+88)=a1
        write(*,*)"system call: ",loopstring(1:l0+88)
        call system(loopstring(1:l0+88))

C read clash.txt (if it exists) and see if it says "CLASH"

        call system('touch clash.txt ; wc -l clash.txt > clash.tmp')
        open(unit=81,file='clash.tmp',status='unknown')
        read(81,*)ijklm
        if(ijklm.ne.0) then
          a6(1:3)='inf'
          write(*,*)'steric clash found in loop conf# ',nmad2-1
          write(*,*)'details follow: '
          write(*,*)
          open(unit=82,file='clash.txt',status='unknown')
678       read(82,'(a66)',end=679)string66
          write(*,'(a66)')string66
          goto 678
679       close(82)
          write(*,*)
        else
          a6(1:3)='   '
          write(*,*)'no steric clash found in loop conf# ',nmad2-1
        endif
        close(81)

C read the energies grepped from loopy's temp1 file 
C note that we used to do the 'ccc' commented line but all lines are now
C commented as we have lost confidence in loopy's ability to determine
C when it's built a hopelessly clashed conformation...

c       call system('rm tempX temp2')
ccc old call system("grep 'energy of vw' temp1 | tail -1 > temp2")
c       call system("grep '^vw: ' temp1 > tempX")

c       if(nmad2.lt.10) then
c         sys_string6(1:31)='tail -X tempX | head -1 > temp2'
c         write(sys_string6( 7: 7),'(i1)')nmad2
c         call system(sys_string6(1:31))
c       else
c         sys_string6(1:32)='tail -XX tempX | head -1 > temp2'
c         write(sys_string6( 7: 8),'(i2)')nmad2
c         call system(sys_string6(1:32))
c       endif

c       open(unit=12,file='temp2',status='unknown')
c       read(12,*,end=616)a1,a2,a3,a4,a5,a6
c       goto 617     ! don't set iloopydone<>0 if temp2 not empty
c616    iloopydone=0 ! =0 if no vw energy was read
c617    continue
c       close(12)

C now we check to see if we have a contiguous loop built...

        loopy_file='loopy_inp_looper_X.pdb'
        write(loopy_file(18:18),'(i1)')nmad2-1

        ntmp_atm=0 ! need to switch format if atoms>100000
        nlink=0
        open(unit=12,file=loopy_file,status='unknown')
741     read(12,'(a80)',end=742)strong
        if(strong(1:4).ne.'ATOM') goto 741
        ntmp_atm=ntmp_atm+1
        if(ntmp_atm.le.99999) then
          read(strong,'(12x,a4,5x,a,i4,4x,3f8.3)')anam,cc,ires,xx,yy,zz
        else
          read(strong,'(13x,a4,5x,a,i4,4x,3f8.3)')anam,cc,ires,xx,yy,zz
        endif

C note that we must have the right residue number and the right chain!

C v3.1  if(anam.eq.' CA '.and.cc.eq.ichn_tot_loops(i).and.
C    &     ires.ge.ibeg_loc_loops(i)+cur_offset(my_frag)-2.and.
C    &     ires.le.iend_loc_loops(i)+cur_offset(my_frag)+2) then
        if(anam.eq.' CA '.and.cc.eq."A".and.
     &     ires.ge.jbeg_loopy-2.and.ires.le.jend_loopy+2) then
          nlink=nlink+1
          xll(nlink)=xx
          yll(nlink)=yy
          zll(nlink)=zz
        endif
        goto 741
742     close(12)
        write(*,*)'Now Checking For Contiguous Loop; nlink = ',nlink

        ifailed=0
        do im=1,nlink-1
          dist2=(xll(im+1)-xll(im))**2+
     &          (yll(im+1)-yll(im))**2+
     &          (zll(im+1)-zll(im))**2
          dist=sqrt(dist2)
          if(dist.gt.4.0) then
            write(*,*)'distance check ',im,dist
            ifailed=1
          endif
        enddo 

C v4.7 BUGFIX - not really a bug but we should also fail this loop if
C nlink = 0 - this indicates that the loop has no CA atom...

        if(nlink.eq.0) ifailed=1

C now make the decision about whether to accept the loop or not

        if(a6(1:3).eq.'inf'.or.iloopydone.eq.0.or.ifailed.eq.1) then
          write(*,*)'FAILURE in loop #',i,' conf# ',nmad2,iloopydone
          if(a6(1:3).eq.'inf') then
c           write(*,*)'infinite loop energy in loop '
            write(*,*)'the loop contains a bad steric clash'
          elseif(iloopydone.eq.0) then
            write(*,*)'loopy timed out and made no pdb file'
          elseif(ifailed.eq.1) then
            write(*,*)'loopy built a ridiculous, non-closed loop'
          endif

C before deciding to increment the loop size we should ask if other loop
C conformations made at the same time are in better shape...

          if(nmade.lt.nloopy_confs) goto 9876

          write(*,*)'looked at all ',nloopy_confs,
     &              ' confs built by loopy now'
          write(*,*)'and none of them are any good :('
          write(*,*)'so we will now attempt to increase loop size'

ccc       if(i.eq.27) stop ! TEMP TEMP TEMP 

C if we failed to build the loop we need to increase its size in the
C hope of building it successfully - here we make use of
C i_in_pdb_fasta_tot which helps us know which of the positions to
C extend the loop in...

          if(i_in_SS_elemnt_tot(ibeg_tot_loops(i)).gt.
     &       i_in_SS_elemnt_tot(iend_tot_loops(i))) then
            iend_loc_loops(i)=iend_loc_loops(i)+1
            iend_tot_loops(i)=iend_tot_loops(i)+1
            write(*,*)'N-ter loop res is in SS element so survives'
            goto 2345
          elseif(i_in_SS_elemnt_tot(ibeg_tot_loops(i)).lt.
     &           i_in_SS_elemnt_tot(iend_tot_loops(i))) then
            ibeg_loc_loops(i)=ibeg_loc_loops(i)-1
            ibeg_tot_loops(i)=ibeg_tot_loops(i)-1
            write(*,*)'C-ter loop res is in SS element so survives'
            goto 2345
          endif

          if(i_in_pdb_fasta_tot(ibeg_tot_loops(i)).eq.
     &       i_in_pdb_fasta_tot(iend_tot_loops(i))) then
            idiff=iend_tot_loops(i)-ibeg_tot_loops(i)+1
            call random_number(u1)
            if(u1.le.0.5) then
              ibeg_loc_loops(i)=ibeg_loc_loops(i)-1
              ibeg_tot_loops(i)=ibeg_tot_loops(i)-1
            else
              iend_loc_loops(i)=iend_loc_loops(i)+1
              iend_tot_loops(i)=iend_tot_loops(i)+1
            endif
          elseif(i_in_pdb_fasta_tot(ibeg_tot_loops(i)).gt.
     &           i_in_pdb_fasta_tot(iend_tot_loops(i))) then
            iend_loc_loops(i)=iend_loc_loops(i)+1
            iend_tot_loops(i)=iend_tot_loops(i)+1
          elseif(i_in_pdb_fasta_tot(ibeg_tot_loops(i)).lt.
     &           i_in_pdb_fasta_tot(iend_tot_loops(i))) then
            ibeg_loc_loops(i)=ibeg_loc_loops(i)-1
            ibeg_tot_loops(i)=ibeg_tot_loops(i)-1
          endif

2345      continue

C we now check whether the extended loop now encroaches onto other 
C loops and/or tail sections - if it does we set it to -999 so that we
C no longer deal with it - I haven't yet fixed it so that loops merge -
C THIS NEEDS TO BE IMPLEMENTED...

C BUGFIX v1.1 we reorder these if statements so that we first check against
C existing tails before considering whether to make a new tail - the
C previous code would fail to merge if the existing tail contained only
C one residue as it would satisfy the first if statement encountered

C make loops that encroach on existing Nter tails part of those tails...
C v2.0 note that we add an additional condition here - it makes sure
C that the loop overlaps but is not contained within the tail - if it's
C totally contained within it then the chances are that it's in a
C structured domain that we'll use the rebuild_tails code for

          do j=1,num_tot_nters
            if(ifrg_tot_loops(i).ne.ifrg_tot_nters(j)) cycle
            if(ibeg_loc_loops(i).le.iend_loc_nters(j).and.
     &         iend_loc_loops(i).ge.iend_loc_nters(j).and. !BUGFIX v2.0
     &         ibeg_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' now merges with Nter tail ',j
              iend_tot_nters(j)=iend_tot_loops(i)
              iend_loc_nters(j)=iend_loc_loops(i)
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          enddo

C make loops that encroach on existing Cter tails part of those tails...
C v2.0 note that we add an additional condition here - it makes sure
C that the loop overlaps but is not contained within the tail - if it's
C totally contained within it then the chances are that it's in a
C structured domain that we'll use the rebuild_tails code for

          do j=1,num_tot_cters
            if(ifrg_tot_loops(i).ne.ifrg_tot_cters(j)) cycle
            if(iend_loc_loops(i).ge.ibeg_loc_cters(j).and.
     &         ibeg_loc_loops(i).le.ibeg_loc_cters(j).and. !BUGFIX v2.0
     &         iend_loc_loops(i).ne.-999) then
              write(*,*)'loop #',i,' now merges with Cter tail ',j
              ibeg_tot_cters(j)=ibeg_tot_loops(i)
              ibeg_loc_cters(j)=ibeg_loc_loops(i)
              ibeg_loc_loops(i)=-999 ! exploit this later
              iend_loc_loops(i)=-999 ! exploit this later
            endif
          enddo

C make loops that are at beg of chain become new Nter tails

          if(ibeg_loc_loops(i).eq.1) then
            num_tot_nters=num_tot_nters+1
            write(*,*)'loop #',i,' is now a new Nter tail ',
     &                num_tot_nters
            ifrg_tot_nters(num_tot_nters)=ifrg_tot_loops(i)
            ibeg_tot_nters(num_tot_nters)=ibeg_tot_loops(i)
            iend_tot_nters(num_tot_nters)=iend_tot_loops(i)

C v2.1 add here

            ichn_tot_nters(num_tot_nters)=ichn_tot_loops(i)

            ibeg_loc_nters(num_tot_nters)=ibeg_loc_loops(i)
            iend_loc_nters(num_tot_nters)=iend_loc_loops(i)
            ibeg_loc_loops(i)=-999 ! exploit this later
            iend_loc_loops(i)=-999 ! exploit this later
          elseif(ibeg_loc_loops(i).lt.1.and.
     &           ibeg_loc_loops(i).ne.-999) then
            write(*,*)'error for beg loop i ',i,ibeg_loc_loops(i)
            stop
          endif

C make loops that are at end of chain become new Cter tails

          if(iend_loc_loops(i).eq.
     &       num_res_in_frag(ifrg_tot_loops(i))) then
            num_tot_cters=num_tot_cters+1
            write(*,*)'loop #',i,' is now a new Cter tail ',
     &                num_tot_cters
            ifrg_tot_cters(num_tot_cters)=ifrg_tot_loops(i)
            ibeg_tot_cters(num_tot_cters)=ibeg_tot_loops(i)
            iend_tot_cters(num_tot_cters)=iend_tot_loops(i)

C v2.1 add here

            ichn_tot_cters(num_tot_cters)=ichn_tot_loops(i)

            ibeg_loc_cters(num_tot_cters)=ibeg_loc_loops(i)
            iend_loc_cters(num_tot_cters)=iend_loc_loops(i)
            ibeg_loc_loops(i)=-999 ! exploit this later
            iend_loc_loops(i)=-999 ! exploit this later
          elseif(iend_loc_loops(i).gt.
     &           num_res_in_frag(ifrg_tot_loops(i))) then
            write(*,*)'error for end loop i ',i,iend_loc_loops(i)
            stop
          endif

          write(*,*)'LOOP ',i,ibeg_tot_loops(i),iend_tot_loops(i)

C now set i_am_original to the newly created ends to be zero - I think
C that this code will still work okay even if we encroach on existing
C tails or if we make a new tail but I need to think about this...
C UPDATE - we should only do this if we haven't already switched this
C loop off by setting it to -999... (also note that in this code we will
C have only changed beg or end by, at most, one residue, which is why we
C only have to look at beg or end to update i_am_original - no need to
C do a loop from beg to end since all intervening res's are already set)

          if(ibeg_loc_loops(i).ne.-999) then
            if(i_am_original(ibeg_tot_loops(i)).ge.1) then
               i_am_original(ibeg_tot_loops(i))=0
            endif
            if(i_am_original(iend_tot_loops(i)).ge.1) then
               i_am_original(iend_tot_loops(i))=0
            endif
          endif

          goto 606

        else
          write(*,*)'SUCCESS in loop #',i,' out of ',num_tot_loops
        endif

C v1.8 - changed mv command to cp here - not sure why we used mv...

ccc     sys_stringy(1:20)='cp loopy_inp_looper_'
ccc     write(sys_stringy(21:21),'(i1)')nmad2-1
ccc     sys_stringy(22:35)='.pdb crap.pdb'
ccc     call system(sys_stringy(1:35))

        loopnam2(1:22)='loopy_inp_looper_X.pdb'
        write(loopnam2(18:18),'(i1)')nmad2-1
        open(unit=11,file=loopnam2,status='unknown') ! input
        open(unit=12,file='crap.pdb',status='unknown') ! output
8561    read(11,'(a80)',end=8562)char80
        write(12,'(a80)')char80
        goto 8561
8562    close(11)
        close(12)

C v2.5 also keep a copy of every loop immediately after making...
C v4.5 - not yet implemented a replacement for system call...

        sys_stringw(1:20)='cp loopy_inp_looper_'
        write(sys_stringw(21:21),'(i1)')nmad2-1
        sys_stringw(22:47)='.pdb loopy_SUCCESS_XXX.pdb'
        write(char3,'(i3)')i
        if(char3(1:1).eq.' ') char3(1:1)='0'
        if(char3(2:2).eq.' ') char3(2:2)='0'
        sys_stringw(41:43)=char3
        call system(sys_stringw(1:47))

C v2.1 - we used to have to correct the loopy_pdb file output for cases
C where the number of atoms exceeded 10,000 but this is almost certainly
C unnecessary now that we only include atoms near the loops - for now,
C we'll keep this in place though...

        sys_stringx(1:len_code)=code_path(1:len_code)
        l0=len_code+1
        sys_stringx(l0:l0+27)='/correct_loopy_pdb_file.exe '
        sys_stringx(l0+28:l0+49)='crap.pdb loopy_out.pdb'
        call system(sys_stringx(1:l0+49))

C at this point we should have loopy_tmp.pdb and loopy_out.pdb - we want
C to add the linker coords found in loopy_out.pdb into loopy_tmp.pdb...
C - first we'll make loopy_nxt.pdb then do a system call 2 loopy_tmp.pdb

        open(unit=11,file='loopy_tmp.pdb',status='unknown')
        open(unit=12,file='loopy_out.pdb',status='unknown')
        open(unit=13,file='loopy_nxt.pdb',status='unknown')

C get the correct values of jbeg and jend at this point

        jbeg=ibeg_loc_loops(i)+cur_offset(my_frag)
        jend=iend_loc_loops(i)+cur_offset(my_frag)

C v3.1 - the repeated setting of jbeg/jend above seems unnecessary so
C the following lines are probably also unnecessary...

        jbeg_loopy=ibeg_loc_loops(i)
        jend_loopy=iend_loc_loops(i)

        iatpoint=0
9007    read(11,'(a66)',end=9008)string66
        read(string66,'(21x,a,i4)')chainID,jkl

C v2.5 here's where we skip intra-loop residues that are present in
C loopy_tmp.pdb but which will be replaced by the newly built loop

        if(chainID.eq.ichn_tot_loops(i).and.
     &     jkl.ge.jbeg.and.
     &     jkl.le.jend) goto 9007

C read until we get to the first atom of the res *after* the loop

        if(chainID.eq.ichn_tot_loops(i).and.jkl.eq.jend+1) then
          if(iatpoint.eq.0) then
9009        read(12,'(a66)',end=9010)strong66
            read(strong66,'(21x,a,i4)')chainID,jkl
C v3.1      if(chainID.eq.ichn_tot_loops(i).and.jkl.ge.jbeg.and.
C    &         jkl.le.jend) write(13,'(a66)')strong66

C v3.1 here's where we convert the chainID and residue numbers back to
C what they should be in the final output

            if(chainID.eq."A".and.jkl.ge.jbeg_loopy.and.
     &         jkl.le.jend_loopy) then
              strong66(22:22)=ichn_tot_loops(i)
              jkl_tmp=jkl+(jbeg-jbeg_loopy)
              write(char4,'(i4)')jkl_tmp
              strong66(23:26)=char4
              write(13,'(a66)')strong66
            endif

            goto 9009
9010        close(12)
            write(13,'(a66)')string66
            iatpoint=1
          else
            write(13,'(a66)')string66
          endif
        else
          write(13,'(a66)')string66
        endif
        goto 9007

9008    close(11)
        close(13)

C v2.5 also keep a copy of every loop immediately after making...

ccc     sys_stringw(1:44)='cp loopy_nxt.pdb loopy_SUCCESS_final_XXX.pdb'
ccc     write(char3,'(i3)')i
ccc     if(char3(1:1).eq.' ') char3(1:1)='0'
ccc     if(char3(2:2).eq.' ') char3(2:2)='0'
ccc     sys_stringw(38:40)=char3
ccc     call system(sys_stringw(1:44))

        loopname(1:27)='loopy_SUCCESS_final_XXX.pdb'
        write(char3,'(i3)')i
        if(char3(1:1).eq.' ') char3(1:1)='0'
        if(char3(2:2).eq.' ') char3(2:2)='0'
        loopname(21:23)=char3

        open(unit=11,file='loopy_nxt.pdb',status='unknown') ! input
        open(unit=12,file=loopname,status='unknown') ! output
8551    read(11,'(a80)',end=8552)char80
        write(12,'(a80)')char80
        goto 8551
8552    close(11)
        close(12)

ccc     call system('mv loopy_nxt.pdb loopy_tmp.pdb')
        open(unit=11,file='loopy_nxt.pdb',status='unknown') ! input
        open(unit=12,file='loopy_tmp.pdb',status='unknown') ! output
8541    read(11,'(a80)',end=8542)char80
        write(12,'(a80)')char80
        goto 8541
8542    close(11)
        close(12)

      enddo ! do i=1,num_tot_loops

C before looking at tails we need to make sure we've got the latest
C version of loopy_inp.pdb made:

ccc   call system('mv loopy_tmp.pdb loopy_inp.pdb')
ccc   write(*,*)'istatus 004 = ',istatus
      open(unit=11,file='loopy_tmp.pdb',status='unknown') ! input
      open(unit=12,file='loopy_inp.pdb',status='unknown') ! output
8531  read(11,'(a80)',end=8532)char80
      write(12,'(a80)')char80
      goto 8531
8532  close(11)
      close(12)


C v3.6 come here if skip_all_loops=.true.

7654  continue 

C at this point, loopy_inp.pdb contains all loops so can now be used as
C a fixed template from which to build the tails - now we need to start
C the hard work of using the dipeptide libraries

C but, before that, if build_no_tails is specified then act accordingly

      if(build_no_tails) then
        write(*,*)
        write(*,*)'truncate_tails = "all" so we will build no tails'
        write(*,*)
        num_tot_nters=0
        num_tot_cters=0
      endif

C if we have any tails to do then read all the dipeptide pdb files:
C v3.4 we will read these files in in the following order

C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_coils.pdb
C OR v4.6
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_norestrict.pdb
C THEN
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Nter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_helix.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Cter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Nter.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_sheet.pdb
C /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Cter.pdb

      if(num_tot_nters.gt.0.or.num_tot_cters.gt.0) then
        write(*,*)
        write(*,*)'time to start reading dipeptide libraries'
        write(*,*)
        allocate(string66_dipep(1:30,1:16000,1:20,1:20))
        allocate(x_dipep(1:30,1:16000,1:20,1:20))
        allocate(y_dipep(1:30,1:16000,1:20,1:20))
        allocate(z_dipep(1:30,1:16000,1:20,1:20))
        allocate(n_temp4(1:10000))

C v4.5 up the max-atoms in each residue from 50 to 10000... this allows
C up to 10,000 HETATM entries to be associated with each residue in the
C tail - this is massive overkill and crude so we may need to reconsider
C this later so that we can allow >10000 with a single residue instead

        allocate(string66_temp4(1:10000,1:10000))
        allocate(d_temp4(1:10000,1:10000)) ! domain# for this res
        allocate(x_temp4(1:10000,1:10000))
        allocate(y_temp4(1:10000,1:10000))
        allocate(z_temp4(1:10000,1:10000))

        do i=1,20
          do j=1,20

            if(i_need_dipep(i,j).eq.0) cycle ! v3.4

            dipepname(1:klm)=dipeptide_path(1:klm)
            dipepname(klm+1:klm+1)='/'
            dipepname(klm+2:klm+2)=dipep(i:i)
            dipepname(klm+3:klm+3)=dipep(j:j)
c           dipepname(klm+4:klm+4)='/'         ! comment out for PISCES
c           dipepname(klm+5:klm+5)=dipep(i:i)  !
c           dipepname(klm+6:klm+6)=dipep(j:j)  !
            if(use_norestrict) then
              dipepname(klm+4:klm+18)='_norestrict.pdb'
            else
              dipepname(klm+4:klm+18)='_both_coils.pdb'
            endif
            write(*,*)'working on ',dipepname(1:klm+18)

C following lines can be used to get num_peps by finding the last entry
C that says MODEL in the _combined.pdb file - comment out once the
C numbers have been saved to a nohup.out - and place the num_peps= lines
C near the top of the program...

c           grepstring(1:10)='tail -30  '
c           grepstring(11:klm+26)=dipepname(1:klm+16)
c         grepstring(klm+27:klm+60)=' | grep MODEL | tail -1 > junk.txt'
c           write(*,*)'calling grepstring ',grepstring(1:klm+60)
c           call system(grepstring(1:klm+60))
c           open(unit=11,file='junk.txt',status='unknown')
c           read(11,'(8x,i6)')num_peps_tmp      
c           num_peps(i,j)=num_peps_tmp
c           close(11)

C the following line is the key one that writes out num_peps in a form
C that can be immediately pasted into the ahemodel code (see above)

c           write(*,'("      num_peps(",i2,",",i2,") = ",i6)')i,j,
c    &                       num_peps(i,j)
 
            open(unit=11,file=dipepname(1:klm+18),status='unknown')

C first, figure out how many atoms in this dipeptide

            natoms=0
            natoms_all=0
            num_N =0
            num_CA=0
            num_C =0
            num_O =0
            read(11,'(1x)')
5991        read(11,'(a66)')string66
            if(string66(1:6).ne.'ENDMDL') then

C use sidechain_mode here to decide how to treat sidechain atoms
C if sidechain_mode=0 then we accept all atoms regardless
C note that natoms_all includes all atoms (duh), while natoms includes
C only those atoms that were accepted according to sidechain_mode...

              natoms_all=natoms_all+1
              if(sidechain_mode.eq.1) then
                if(string66(13:16).ne.' N  '.and. 
     &             string66(13:16).ne.' CA '.and.
     &             string66(13:16).ne.' C  '.and.
     &             string66(13:16).ne.' O  '.and.
     &             string66(13:16).ne.' CB ') goto 5991
              elseif(sidechain_mode.eq.2) then
                if(string66(13:16).ne.' N  '.and. 
     &             string66(13:16).ne.' CA '.and.
     &             string66(13:16).ne.' C  '.and.
     &             string66(13:16).ne.' O  ') goto 5991
              endif

              natoms=natoms+1
              if(natoms.eq.1)string66_first=string66
              if(string66(18:26).ne.string66_first(18:26)) then
                string66_first=string66 ! ensures we don't keep doing
                dipep_natom1(i,j)=natoms-1
              endif
              if(string66(13:16).eq.' N  ') then
                num_N=num_N+1
                if(num_N.eq.1) then
                  dipep_N1(i,j)=natoms
                else
                  dipep_N2(i,j)=natoms
                endif
              elseif(string66(13:16).eq.' CA ') then
                num_CA=num_CA+1
                if(num_CA.eq.1) then
                  dipep_CA1(i,j)=natoms
                else
                  dipep_CA2(i,j)=natoms
                endif
              elseif(string66(13:16).eq.' C  ') then
                num_C=num_C+1
                if(num_C.eq.1) then
                  dipep_C1(i,j)=natoms
                else
                  dipep_C2(i,j)=natoms
                endif
              elseif(string66(13:16).eq.' O  ') then
                num_O=num_O+1
                if(num_O.eq.1) then
                  dipep_O1(i,j)=natoms
                else
                  dipep_O2(i,j)=natoms
                endif
              endif
            else
              goto 5992
            endif
            goto 5991
5992        rewind(11)
            dipep_natoms(i,j)=natoms
            dipep_natom2(i,j)=natoms-dipep_natom1(i,j)
c           write(*,*)'dipep_natoms = ',dipep_natoms(i,j),
c    &      dipep_natom1(i,j),dipep_natom2(i,j),num_peps(i,j),
c    &      dipep_N1(i,j),dipep_N2(i,j),dipep_CA1(i,j),dipep_CA2(i,j),
c    &      dipep_C1(i,j),dipep_C2(i,j),dipep_O1(i,j),dipep_O2(i,j)

C now read all models...

            dipep_both_coils_beg(i,j)=1
            nmods=0
            k=0 ! this will keep track of total confs for this dipep
            do kf=1,100000
              read(11,'(1x)',end=47801) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ') 
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47801       close(11)
            dipep_both_coils_end(i,j)=dipep_both_coils_beg(i,j)+nmods-1

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Nter.pdb

C v3.9 here we have an option to skip the reading of helix dipeptide
C files - we only read them if we both: (1) have at least one ss2 file
C and (2) if the helix_cut variable is <= 1...

C NOTE THAT I HAVE NOT TESTED THE SITUATION WHERE ONLY ONE OF HELIX_CUT
C OR SHEET_CUT IS SET >1 - WE SHOULD SET BOTH OR NONE TO >1 FOR NOW...

            if((helix_cut.gt.1.0.and.sheet_cut.le.1.0).or.
     &         (helix_cut.le.1.0.and.sheet_cut.gt.1.0)) then
              write(*,*)
              write(*,*)'if one or other of helix_cut or sheet_cut '
              write(*,*)'is set to >1 then *both* must be set...'
              write(*,*)'since the current settings have not been '
              write(*,*)'tested explicitly, this is set to a fatal '
              write(*,*)'error for now - so either set both to >1 '
              write(*,*)'or set both to <=1...'
              write(*,*)
              stop
            endif

            if(have_ss2_file.and.helix_cut.le.1.0) then 

            dipepname(klm+4:klm+18)='_helix_Nter.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
            dipep_helix_Nter_beg(i,j)=dipep_both_coils_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47802) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ') 
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47802       close(11)
            dipep_helix_Nter_end(i,j)=dipep_helix_Nter_beg(i,j)+nmods-1

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_helix.pdb

            dipepname(klm+4:klm+18)='_both_helix.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
            dipep_both_helix_beg(i,j)=dipep_helix_Nter_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47803) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ')
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47803       close(11)
            dipep_both_helix_end(i,j)=dipep_both_helix_beg(i,j)+nmods-1

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_helix_Cter.pdb

            dipepname(klm+4:klm+18)='_helix_Cter.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
            dipep_helix_Cter_beg(i,j)=dipep_both_helix_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47804) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ')
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47804       close(11)
            dipep_helix_Cter_end(i,j)=dipep_helix_Cter_beg(i,j)+nmods-1
            endif ! v3.9 skip if no helix wanted

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Nter.pdb

C v3.9 here we have an option to skip the reading of sheet dipeptide
C files - we only read them if we both: (1) have at least one ss2 file
C and (2) if the sheet_cut variable is <= 1...

            if(have_ss2_file.and.sheet_cut.le.1.0) then 

            dipepname(klm+4:klm+18)='_sheet_Nter.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
C 2023 BUG! dipep_sheet_Nter_beg(i,j)=dipep_both_coils_end(i,j)+1
            dipep_sheet_Nter_beg(i,j)=dipep_helix_Cter_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47805) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ') 
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47805       close(11)
            dipep_sheet_Nter_end(i,j)=dipep_sheet_Nter_beg(i,j)+nmods-1

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_both_sheet.pdb

            dipepname(klm+4:klm+18)='_both_sheet.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
            dipep_both_sheet_beg(i,j)=dipep_sheet_Nter_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47806) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ')
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47806       close(11)
            dipep_both_sheet_end(i,j)=dipep_both_sheet_beg(i,j)+nmods-1

C now do /home/aelcock/2017_ECOLI_CELL/PISCES/WF_sheet_Cter.pdb

            dipepname(klm+4:klm+18)='_sheet_Cter.pdb'
            write(*,*)'working on ',dipepname(1:klm+18)
            open(unit=11,file=dipepname(1:klm+18),status='unknown')
            dipep_sheet_Cter_beg(i,j)=dipep_both_sheet_end(i,j)+1
            nmods=0
            do kf=1,100000
              read(11,'(1x)',end=47807) ! skip the 'MODEL' line
              nmods=nmods+1
              k=k+1
              l=0
              do n=1,natoms_all
                read(11,'(a66)')string66

                if(sidechain_mode.eq.1) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  '.and.
     &               string66(13:16).ne.' CB ') cycle
                elseif(sidechain_mode.eq.2) then
                  if(string66(13:16).ne.' N  '.and. 
     &               string66(13:16).ne.' CA '.and.
     &               string66(13:16).ne.' C  '.and.
     &               string66(13:16).ne.' O  ') cycle
                endif
                l=l+1

                string66_dipep(l,k,i,j)=string66
                read(string66,'(30x,3f8.3)')x,y,z
                x_dipep(l,k,i,j)=x
                y_dipep(l,k,i,j)=y
                z_dipep(l,k,i,j)=z
                if(l.eq.dipep_N1(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_N2(i,j).and.string66(13:16).ne.' N  '.or.
     &             l.eq.dipep_CA1(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_CA2(i,j).and.string66(13:16).ne.' CA '.or.
     &             l.eq.dipep_C1(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_C2(i,j).and.string66(13:16).ne.' C  '.or.
     &             l.eq.dipep_O1(i,j).and.string66(13:16).ne.' O  '.or.
     &             l.eq.dipep_O2(i,j).and.string66(13:16).ne.' O  ')
     &          then
                  write(*,*)'oops ',dipep(i:i),' ',dipep(j:j),' ',
     &                      'k = ',k,' l = ',l,' not N,CA,C,O ',
     &                      string66(13:16)
                endif
              enddo
              read(11,'(1x)') ! skip the 'ENDMDL' line
            enddo
47807       close(11)
            dipep_sheet_Cter_end(i,j)=dipep_sheet_Cter_beg(i,j)+nmods-1
            endif ! v3.9 skip if no sheet wanted

            num_peps(i,j)=k

c           write(*,47815)i,j,dipep_both_coils_beg(i,j),
c    &                        dipep_both_coils_end(i,j),
c    &                        dipep_helix_Nter_beg(i,j),
c    &                        dipep_helix_Nter_end(i,j),
c    &                        dipep_both_helix_beg(i,j),
c    &                        dipep_both_helix_end(i,j),
c    &                        dipep_helix_Cter_beg(i,j),
c    &                        dipep_helix_Cter_end(i,j)
47815       format('checkage ',11i6)

          enddo
        enddo
        write(*,*)
        write(*,*)'done with reading dipeptide libraries'
        write(*,*)
      endif

C do a fake allocation for temp_grid and real_grid

      nx=10
      ny=10
      nz=10
      np=1
      allocate(temp_grid(1:nx,1:ny,1:nz))
      allocate(cchn_grid(1:np,1:nx,1:ny,1:nz))
      allocate(rchn_grid(1:np,1:nx,1:ny,1:nz))
      allocate(ihet_grid(1:np,1:nx,1:ny,1:nz)) ! = 1 if HETATM
      allocate(xatm_grid(1:np,1:nx,1:ny,1:nz))
      allocate(yatm_grid(1:np,1:nx,1:ny,1:nz))
      allocate(zatm_grid(1:np,1:nx,1:ny,1:nz))

C before going on we store the current definitions of tails so that we
C can restore them when building additional models...
 
ccc   call system('cp loopy_inp.pdb loopy_inp.pdb.original')
ccc   write(*,*)'istatus 003 = ',istatus
      open(unit=11,file='loopy_inp.pdb',status='unknown') ! input
      open(unit=12,file='loopy_inp.pdb.original',status='unknown') ! output
8521  read(11,'(a80)',end=8522)char80
      write(12,'(a80)')char80
      goto 8521
8522  close(11)
      close(12)

      i_am_original_original=i_am_original
      inot_made_yet_original=inot_made_yet

      do i=1,num_tot_nters
        ibeg_loc_nters_original(i)=ibeg_loc_nters(i)
        iend_loc_nters_original(i)=iend_loc_nters(i)
        ibeg_tot_nters_original(i)=ibeg_tot_nters(i)
        iend_tot_nters_original(i)=iend_tot_nters(i)
        ifrg_tot_nters_original(i)=ifrg_tot_nters(i)
        ichn_tot_nters_original(i)=ichn_tot_nters(i)
      enddo

      do i=1,num_tot_cters
        ibeg_loc_cters_original(i)=ibeg_loc_cters(i)
        iend_loc_cters_original(i)=iend_loc_cters(i)
        ibeg_tot_cters_original(i)=ibeg_tot_cters(i)
        iend_tot_cters_original(i)=iend_tot_cters(i)
        ifrg_tot_cters_original(i)=ifrg_tot_cters(i)
        ichn_tot_cters_original(i)=ichn_tot_cters(i)
      enddo

C come back here to start building multiple AHEMODEL structures
C NOT SURE YET IF THIS WORKS WITH REBUILD TAILS - WORTH TRYING

      nmodels=0

C v4.7 let's start timing...

      write(*,*)
      write(*,*)'beginning timing the whole process now'
      write(*,*)
      call system_clock(count=ktime1,count_rate=ru)

1234  continue

C restore the original definitions...

ccc   call system('cp loopy_inp.pdb.original loopy_inp.pdb')
ccc   write(*,*)'istatus 002 = ',istatus
      open(unit=11,file='loopy_inp.pdb.original',status='unknown') ! input
      open(unit=12,file='loopy_inp.pdb',status='unknown') ! output
8511  read(11,'(a80)',end=8512)char80
      write(12,'(a80)')char80
      goto 8511
8512  close(11)
      close(12)

      call system('echo "HI THERE 001"')

      i_am_original=i_am_original_original
      inot_made_yet=inot_made_yet_original

      do i=1,num_tot_nters
        ibeg_loc_nters(i)=ibeg_loc_nters_original(i)
        iend_loc_nters(i)=iend_loc_nters_original(i)
        ibeg_tot_nters(i)=ibeg_tot_nters_original(i)
        iend_tot_nters(i)=iend_tot_nters_original(i)
        ifrg_tot_nters(i)=ifrg_tot_nters_original(i)
        ichn_tot_nters(i)=ichn_tot_nters_original(i)
      enddo

      do i=1,num_tot_cters
        ibeg_loc_cters(i)=ibeg_loc_cters_original(i)
        iend_loc_cters(i)=iend_loc_cters_original(i)
        ibeg_tot_cters(i)=ibeg_tot_cters_original(i)
        iend_tot_cters(i)=iend_tot_cters_original(i)
        ifrg_tot_cters(i)=ifrg_tot_cters_original(i)
        ichn_tot_cters(i)=ichn_tot_cters_original(i)
      enddo

C -------------------------------------------
C START BUILDING THE N-TERMINAL TAILS HERE...
C -------------------------------------------

      inot_made_yet=1 ! not sure if this is necessary...

C first, make a combined list of all tails - Nter is type 1; Cter is
C type 2...

      num_tot_tails=0
      do i=1,num_tot_nters
        num_tot_tails=num_tot_tails+1
        typ_tot_tails(num_tot_tails)=1
      enddo
      do i=1,num_tot_cters
        num_tot_tails=num_tot_tails+1
        typ_tot_tails(num_tot_tails)=2
      enddo
      
C start what is in effect a loop over all of the tails that we need to
C build and select them in a random order - this should ensure no bias
C in the ocnstruction of tails

      imadtail(1:num_tot_tails)=0
      jmadtail(1:num_tot_tails)=0
      do nmad=1,num_tot_tails
9765    call random_number(u1)
        itail=int(u1*num_tot_tails)+1
ccc     itail=nmad ! use this to do Nters in order, then Cters in order
        if(imadtail(itail).eq.0) then
          jmadtail(nmad)=itail
          imadtail(itail)=1
        else
          goto 9765
        endif
        write(*,*)'check randomized tail selection ',nmad,
     &            jmadtail(nmad),typ_tot_tails(jmadtail(nmad)),u1
      enddo

C v2.6: start a loop over all tails

      do iii=1,num_tot_tails

C find where in the list of all tails the current one lies

        itmp=jmadtail(iii)

C now, depending on the type of tail (N-ter vs. C-ter) find its number
C in that list of tails - we'll use "i" here since that was used as the
C do-loop counter in older code versions that looped separately over the
C N-terminal tails and C-terminal tails - should simplify code rewriting

        if(typ_tot_tails(itmp).eq.1) then
          i=itmp
        elseif(typ_tot_tails(itmp).eq.2) then
          i=itmp-num_tot_nters
        endif

C v2.1 keep the fragment# of this tail in a temporary variable

        if(typ_tot_tails(itmp).eq.1) then
          my_frag=ifrg_tot_nters(i)
          write(*,*)
          write(*,*)'working on N-ter tail ',i
          write(*,*)'starts at res# ',ibeg_tot_nters(i)
          write(*,*)'ends   at res# ',iend_tot_nters(i)
          write(*,*)'is from frag#  ',ifrg_tot_nters(i)
          write(*,*)'and is chainID ',ichn_tot_nters(i)
          write(*,*)'and is cur_chainID ',cur_chainID(my_frag)
          write(*,*)

C v4.4 set boundary value once for the relaxed residues - if j>jrelax
C it'll be given relaxed clash criteria

          jrelax=iend_tot_nters(i)-nrelax

        elseif(typ_tot_tails(itmp).eq.2) then
          my_frag=ifrg_tot_cters(i)
          write(*,*)
          write(*,*)'working on C-ter tail ',i
          write(*,*)'starts at res# ',ibeg_tot_cters(i)
          write(*,*)'ends   at res# ',iend_tot_cters(i)
          write(*,*)'is from frag#  ',ifrg_tot_cters(i)
          write(*,*)'and is chainID ',ichn_tot_cters(i)
          write(*,*)'and is cur_chainID ',cur_chainID(my_frag)
          write(*,*)

C v4.4 set boundary value once for the relaxed residues - if j<jrelax
C it'll be given relaxed clash criteria

          jrelax=ibeg_tot_cters(i)+nrelax

        endif

        isentback=0  ! use this to ensure we don't immediately extend
                     ! the chain - when we need to extend isentback=1
404     continue     ! come here if we need to extend this tail
        imin_chn=999999  ! min res# found for N-ter builds
        imax_chn=-1      ! max res# found for C-ter builds
        nfails_chn=0 ! #failures for this entire chain
        itochk=0     ! v2.1 use to check intra-tail

C remove any tail pdb file that might still be kicking around

ccc     call system('rm current_tail.pdb')

C if isentback=1 then we must have been sent back to extend the chain

        if(isentback.eq.1) then
          if(typ_tot_tails(itmp).eq.1) then

            if(iend_tot_nters(i).eq.
     &         my_last_res(ifrg_tot_nters(i))) then
              write(*,*)
              write(*,*)'cannot extend Nter chain past last residue'
              write(*,*)'this is fatal - quitting :('
              write(*,*)
              stop
            endif

            iend_tot_nters(i)=iend_tot_nters(i)+1
            iend_loc_nters(i)=iend_loc_nters(i)+1 ! was mssing b4 v2.7
            if(i_am_original(iend_tot_nters(i)).ge.1) then
               i_am_original(iend_tot_nters(i))=0
            endif

C we should also make sure that we don't use a structured version of
C this residue as this was probably the problem to start with...

            write(*,*)
            write(*,*)'setting "inot_made_yet" to 1 for res# ',
     &                 iend_tot_nters(i)
            write(*,*)'new ibeg & iend values ',ibeg_tot_nters(i),
     &                                          iend_tot_nters(i)
            write(*,*)
            inot_made_yet(iend_tot_nters(i))=1

C v4.2 I think we should also make sure that the preceding residue is
C flagged as not made also...

            write(*,*)'v4.2 ALSO setting "inot_made_yet" to 1 for res#',
     &                 iend_tot_nters(i)-1
            write(*,*)
            inot_made_yet(iend_tot_nters(i)-1)=1

          elseif(typ_tot_tails(itmp).eq.2) then

            if(ibeg_tot_cters(i).eq.
     &         my_frst_res(ifrg_tot_cters(i))) then
              write(*,*)
              write(*,*)'cannot extend Cter chain past frst residue'
              write(*,*)'this is fatal - quitting :('
              write(*,*)
              stop
            endif

            ibeg_tot_cters(i)=ibeg_tot_cters(i)-1
            ibeg_loc_cters(i)=ibeg_loc_cters(i)-1 ! was mssing b4 v2.7
            if(i_am_original(ibeg_tot_cters(i)).ge.1) then
               i_am_original(ibeg_tot_cters(i))=0
            endif

C we should also make sure that we don't use a structured version of
C this residue as this was probably the problem to start with...

            write(*,*)
            write(*,*)'setting "inot_made_yet" to 1 for res# ',
     &                 ibeg_tot_cters(i)
            write(*,*)'new ibeg & iend values ',ibeg_tot_cters(i),
     &                                          iend_tot_cters(i)
            write(*,*)
            inot_made_yet(ibeg_tot_cters(i))=1

C v4.2 I think we should also make sure that the succeeding residue is
C flagged as not made also...

            write(*,*)'v4.2 ALSO setting "inot_made_yet" to 1 for res#',
     &                 ibeg_tot_cters(i)+1
            write(*,*)
            inot_made_yet(ibeg_tot_cters(i)+1)=1

          endif
        endif

C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

C v2.1 at this point we can put all nearby atoms onto a grid...
C note that we come back here because this is *after* 404 - i.e. the
C line that we go to when we extend a chain...

        deallocate(temp_grid)
        deallocate(cchn_grid)
        deallocate(rchn_grid)
        deallocate(ihet_grid)
        deallocate(xatm_grid)
        deallocate(yatm_grid)
        deallocate(zatm_grid)

        call system('echo "HI THERE 101"')

C v2.1 set residue numbers correctly here - anything that depends on a
C correctly numbered residue from here on out will use "jx" not "j"

        if(typ_tot_tails(itmp).eq.1) then
          j=iend_tot_nters(i)+1
          jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1
          write(*,*)'TEMP LOOKING FOR JX (1ST RES AFTER TAIL) = ',jx
          write(*,*)'TEMP LOOKING FOR CH = ',ichn_tot_nters(i)
        elseif(typ_tot_tails(itmp).eq.2) then
          j=ibeg_tot_cters(i)-1
          jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1
          write(*,*)'TEMP LOOKING FOR JX (1ST RES BEFOR TAIL) = ',jx
          write(*,*)'TEMP LOOKING FOR CH = ',ichn_tot_cters(i)
        endif

        open(unit=11,file='loopy_inp.pdb',status='unknown')
9011    read(11,'(a66)',end=9012)string66
        read(string66,'(21x,a,i4)')chainID,jkl
        if(typ_tot_tails(itmp).eq.1) then
          if(chainID.eq.ichn_tot_nters(i).and.jkl.eq.jx) then
            if(string66(13:16).eq.' CA ') then
              read(string66,'(30x,3f8.3)')xt,yt,zt
              goto 9012
            endif
          endif
        elseif(typ_tot_tails(itmp).eq.2) then
          if(chainID.eq.ichn_tot_cters(i).and.jkl.eq.jx) then
            if(string66(13:16).eq.' CA ') then
              read(string66,'(30x,3f8.3)')xt,yt,zt
              goto 9012
            endif
          endif
        endif
        goto 9011
9012    rewind(11)



C 2023_v1.1 - here is where we introduce a much more reasonable setup
C for a clash grid - we'll just read all pre-placed atoms (including any
C loops/tails that were added on in previous rounds) - and we'll make
C the clash grid fit snugly around these atoms...

        xmin=+999999.9
        ymin=+999999.9
        zmin=+999999.9
        xmax=-999999.9
        ymax=-999999.9
        zmax=-999999.9

8013    read(11,'(a66)',end=8014)string66

C v4.2 read in occupancy - if -99.99 then it's a HETATM

        read(string66,'(21x,a,i4,4x,3f8.3,f6.2)')chainID,jkl,xt,yt,zt,ot

C TEMP TEMP TEMP following lines guarantee we keep these atoms

        if(ot.eq.99.99) goto 8113
        if(jkl.le.0) goto 8113

C if the atom is already in loopy_inp.pdb but will be replaced by the 
C new tail (e.g. if we're cutting into structure), then skip it too
C note that we need to convert the res# to its equivalent - this is
C awkward to do but doable - we loop over all fragments, find one with
C the same chainID as the current atom, convert to an index
C (chainA->1,B->2 etc), then use that plus the res# read from the pdb
C (which will max out at 9999) to get the true total res# (jxtmp) using
C the 2D converting array "chainlocal_2_totres"...

        jxtmp=0
        do nj=1,chainIDmax
          if(chainID.eq.blah(nj:nj)) then
            nxxx=nj
            exit
          endif
        enddo
        jxtmp=chainlocal_2_totres(nxxx,jkl)

        if(typ_tot_tails(itmp).eq.1) then
          if(chainID.eq.ichn_tot_nters(i).and.
     &         jxtmp.ge.ibeg_tot_nters(i).and.
     &         jxtmp.le.iend_tot_nters(i)) then
c           write(*,*)'not putting res# ',jkl,' on the grid'
c           write(*,*)'in loopy_inp.pdb this appears as res# ',jxtmp,
c    &                ' and chainID ',chainID
            goto 8013
          endif
        elseif(typ_tot_tails(itmp).eq.2) then
          if(chainID.eq.ichn_tot_cters(i).and.
     &         jxtmp.ge.ibeg_tot_cters(i).and.
     &         jxtmp.le.iend_tot_cters(i)) then
c           write(*,*)'not putting res# ',jkl,' on the grid'
c           write(*,*)'in loopy_inp.pdb this appears as res# ',jxtmp,
c    &                ' and chainID ',chainID
            goto 8013
          endif
        endif

C if we get here then let's update the min/max coordinates

8113    continue

        xmin=min(xmin,xt)
        ymin=min(ymin,yt)
        zmin=min(zmin,zt)
        xmax=max(xmax,xt)
        ymax=max(ymax,yt)
        zmax=max(zmax,zt)

        goto 8013

8014    rewind(11)

        write(*,*)
        write(*,*)'xmin/xmax ',xmin,xmax
        write(*,*)'ymin/ymax ',ymin,ymax
        write(*,*)'zmin/zmax ',zmin,zmax
        write(*,*)

        xmin=xmin-5.0
        ymin=ymin-5.0
        zmin=zmin-5.0
        xmax=xmax+5.0
        ymax=ymax+5.0
        zmax=zmax+5.0

C set cell size to be max of any clash distance...

        xcll=max(clash_dist1,clash_dist2)
        ycll=max(clash_dist1,clash_dist2)
        zcll=max(clash_dist1,clash_dist2)

C ...but don't let it get smaller than 2A...

        xcll=max(xcll,2.0)
        ycll=max(ycll,2.0)
        zcll=max(zcll,2.0)

        nx=int((xmax-xmin)/xcll)+1
        ny=int((ymax-ymin)/ycll)+1
        nz=int((zmax-zmin)/zcll)+1
        xinv=1.0/xcll
        yinv=1.0/ycll
        zinv=1.0/zcll

        write(*,*)'construct nonbonded grid w/cen @ ',xt,yt,zt
        write(*,*)'REALLY must check dimensions 1 ',nx,ny,nz,xcll

        allocate(temp_grid(1:nx,1:ny,1:nz))
        temp_grid=0
 
        ninbox_max=0
 
9013    read(11,'(a66)',end=9014)string66
ccc     write(*,*)string66

C v4.2 read in occupancy - if -99.99 then it's a HETATM

ccc     read(string66,'(21x,a,i4,4x,3f8.3)')chainID,jkl,xt,yt,zt
        read(string66,'(21x,a,i4,4x,3f8.3,f6.2)')chainID,jkl,xt,yt,zt,ot

C if an atom is off the grid skip it

        if(xt.lt.xmin.or.xt.gt.xmax.or.
     &     yt.lt.ymin.or.yt.gt.ymax.or.
     &     zt.lt.zmin.or.zt.gt.zmax) goto 9013

C if it's a negative res# then automatically put on the grid...
C v4.2 change this to if ot=-99.99
C v4.5REDUX - TEMP TEMP TEMP - should the following say -99.99???
C I need to remind myself of whether we *intend* to keep HETATMs for the
C nonbonded grid or not...

        if(ot.eq.99.99) goto 9113
        if(jkl.le.0) goto 9113

C if the atom is already in loopy_inp.pdb but will be replaced by the 
C new tail (e.g. if we're cutting into structure), then skip it too
C note that we need to convert the res# to its equivalent - this is
C awkward to do but doable - we loop over all fragments, find one with
C the same chainID as the current atom, convert to an index
C (chainA->1,B->2 etc), then use that plus the res# read from the pdb
C (which will max out at 9999) to get the true total res# (jxtmp) using
C the 2D converting array "chainlocal_2_totres"...

        jxtmp=0
        do nj=1,chainIDmax
          if(chainID.eq.blah(nj:nj)) then
            nxxx=nj
            exit
          endif
        enddo
        jxtmp=chainlocal_2_totres(nxxx,jkl)

        if(typ_tot_tails(itmp).eq.1) then
          if(chainID.eq.ichn_tot_nters(i).and.
     &         jxtmp.ge.ibeg_tot_nters(i).and.
     &         jxtmp.le.iend_tot_nters(i)) then
c           write(*,*)'not putting res# ',jkl,' on the grid'
c           write(*,*)'in loopy_inp.pdb this appears as res# ',jxtmp,
c    &                ' and chainID ',chainID
            goto 9013
          endif
        elseif(typ_tot_tails(itmp).eq.2) then
          if(chainID.eq.ichn_tot_cters(i).and.
     &         jxtmp.ge.ibeg_tot_cters(i).and.
     &         jxtmp.le.iend_tot_cters(i)) then
c           write(*,*)'not putting res# ',jkl,' on the grid'
c           write(*,*)'in loopy_inp.pdb this appears as res# ',jxtmp,
c    &                ' and chainID ',chainID
            goto 9013
          endif
        endif

9113    continue

        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              ninbox_max=max(ninbox_max,temp_grid(i1,j1,k1))
            enddo
          enddo
        enddo
        goto 9013
9014    rewind(11)

C v4.5REDUX I believe the following lines are now bullshit so remove?
C ---------------------------------------------------------------------
C TEMP TEMP TEMP - eventually we need to decide how to treat HETATM
C entries that are associated with structured domains that might be
C being rebuilt - somehow their coordinates need to follow their
C "parent" domain - in the meantime we need to remember that we should
C also probably be adding HETATM entries to the nonbonded grid - so we
C need something like:

c       do m1=1,num_het_atms
c         xt=xh(m1)
c         yt=yh(m1)
c         zt=zh(m1)

C put me on the grid etc...

c       enddo
C ---------------------------------------------------------------------

        write(*,*)'ninbox_max = ',ninbox_max
        write(*,*)'allocating memory for grid'

        allocate(cchn_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)
        allocate(rchn_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)
        allocate(ihet_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)
        allocate(xatm_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)
        allocate(yatm_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)
        allocate(zatm_grid(1:ninbox_max,1:nx,1:ny,1:nz),stat=istat)

        temp_grid=0

C v4.2 set counters for ATOM-ATOM clashers and ATOM-HETATM clashers

        naa=0
        nha=0

9015    read(11,'(a66)',end=9016)string66

C v4.2 read in occupancy - if -99.99 then it's a HETATM

ccc     read(string66,'(21x,a,i4,4x,3f8.3)')chainID,jkl,xt,yt,zt
        read(string66,'(21x,a,i4,4x,3f8.3,f6.2)')chainID,jkl,xt,yt,zt,ot

        if(xt.lt.xmin.or.xt.gt.xmax.or.
     &     yt.lt.ymin.or.yt.gt.ymax.or.
     &     zt.lt.zmin.or.zt.gt.zmax) goto 9015

C if it's a negative res# then automatically put on the grid...
C v4.2 change this to if ot=-99.99
C v4.5REDUX - TEMP TEMP TEMP - should the following say -99.99???
C I need to remind myself of whether we *intend* to keep HETATMs for the
C nonbonded grid or not...

        if(ot.eq.99.99) goto 9116
        if(jkl.le.0) goto 9116

C if the atom is already in loopy_inp.pdb but will be replaced by the 
C new tail (e.g. if we're cutting into structure), then skip it too
C note that we need to convert the res# to its equivalent as above

        jxtmp=0
        do nj=1,chainIDmax
          if(chainID.eq.blah(nj:nj)) then
            nxxx=nj
            exit
          endif
        enddo
        jxtmp=chainlocal_2_totres(nxxx,jkl)

        if(typ_tot_tails(itmp).eq.1) then
          if(chainID.eq.ichn_tot_nters(i).and.
     &         jxtmp.ge.ibeg_tot_nters(i).and.
     &         jxtmp.le.iend_tot_nters(i)) then
            goto 9015
          endif
        elseif(typ_tot_tails(itmp).eq.2) then
          if(chainID.eq.ichn_tot_cters(i).and.
     &         jxtmp.ge.ibeg_tot_cters(i).and.
     &         jxtmp.le.iend_tot_cters(i)) then
            goto 9015
          endif
        endif

9116    continue

        i1=int((xt-xmin)*xinv)+1
        j1=int((yt-ymin)*yinv)+1
        k1=int((zt-zmin)*zinv)+1
        ilo=i1-1
        jlo=j1-1
        klo=k1-1
        ihi=i1+1
        jhi=j1+1
        khi=k1+1
        if(ilo.lt.1) ilo=1
        if(jlo.lt.1) jlo=1
        if(klo.lt.1) klo=1
        if(ihi.gt.nx) ihi=nx
        if(jhi.gt.ny) jhi=ny
        if(khi.gt.nz) khi=nz
        ido=0 ! use to only increment naa and nha once per atom
        do k1=klo,khi
          do j1=jlo,jhi
            do i1=ilo,ihi
              temp_grid(i1,j1,k1)=temp_grid(i1,j1,k1)+1
              cchn_grid(temp_grid(i1,j1,k1),i1,j1,k1)=chainID
              if(ot.eq.-99.99) then
                ihet_grid(temp_grid(i1,j1,k1),i1,j1,k1)=1
                if(ido.eq.0) nha=nha+1
                ido=1
              else
                ihet_grid(temp_grid(i1,j1,k1),i1,j1,k1)=0
                if(ido.eq.0) naa=naa+1
                ido=1
              endif
              rchn_grid(temp_grid(i1,j1,k1),i1,j1,k1)=jkl
              xatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=xt
              yatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=yt
              zatm_grid(temp_grid(i1,j1,k1),i1,j1,k1)=zt
            enddo
          enddo
        enddo
        goto 9015
9016    close(11)
        write(*,*)'done putting atoms on a grid'
        write(*,*)'#ATOM   entries on grid = ',naa
        write(*,*)'#HETATM entries on grid = ',nha

C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

C v4.3/v4.4 store some stuff

        ngiveups=0    

5432    continue ! come here if giving up

        itrybacktrack=0
        nbackmoves=0
        nbackresets=0
        itochk=0 ! I think we need this if we are sent back to 5432

        if(typ_tot_tails(itmp).eq.1) then
          jstart    =iend_tot_nters(i)
          jbacktrack=iend_tot_nters(i)
        elseif(typ_tot_tails(itmp).eq.2) then
          jstart    =ibeg_tot_cters(i)
          jbacktrack=ibeg_tot_cters(i)
        endif

505     continue ! come here if we're restarting this tail

C remove any tail pdb file that might still be kicking around
C also need to reset itochk to zero...

ccc     call system('rm current_tail.pdb')

C v4.3/v4.4
 
        if(itrybacktrack.eq.1) then

C do N-terminal tails here

          if(typ_tot_tails(itmp).eq.1) then

C v4.4 check if we've failed nbacktries times - if we have then we'll
C try to reset jbacktrack further back up the chain

            nbackmoves=nbackmoves+1
            if(nbackmoves.ge.nbacktries) then

C if we've repeatedly restarted it may be time to start from scratch

              nbackresets=nbackresets+1
              if(nbackresets.ge.nbacktries) then

C first check if we need to delete structure from the tail - note that
C we only do this if the logical "iworkoutwards" is invoked - if it is
C not then we do nothing and just restart the tail

                ngiveups=ngiveups+1
                if(ngiveups.ge.nbacktries.and.iworkoutwards) then
                  kwant=0
                  do k=iend_tot_nters(i),ibeg_tot_nters(i),-1
                    if(inot_made_yet(k).eq.0) then
                      kwant=k
                      exit
                    endif
                  enddo
                  if(kwant.ne.0) then
                    write(*,*)'giving up^2; remove struc frm res#',kwant
                    inot_made_yet(kwant)=1
                  endif
                  ngiveups=0
                  goto 5432
                else
                  write(*,*)'giving up; start from scratch # ',ngiveups
                  goto 5432
                endif

              endif

C v4.4 count up to find another backtrack position: jbacktrack - note
C that we only count flexible residues in the count...

              jbacktrack_old=jbacktrack
              imin_chn=jbacktrack_old  ! reset min res# found for N-ter builds
              l=0
              kwant=0
              do k=jbacktrack+1,iend_tot_nters(i)
                if(inot_made_yet(k).eq.1) then
                  l=l+1
                  if(l.eq.nbacktrack) then
                    kwant=k
                    exit
                  endif
                endif
              enddo
              if(kwant.ne.0) then
                jbacktrack=kwant
              else
                jbacktrack=iend_tot_nters(i)
              endif
              nbackmoves=0
              if(jbacktrack.ne.jbacktrack_old) then
              write(*,*)'NEW backtrack position for model# ',nmodels+1,
     &                  ' jbacktrack = ',jbacktrack
              endif
            endif

            jstart=jbacktrack

C now keep only those atoms for which j>jstart
C v4.5 this will need checking to see if it works if we have HETATM
C entries - I think that as written it won't - TEMP TEMP TEMP - so for
C now we should only try HETATM-association with nbacktrack=0
 
            itochk_tmp=0
            do ito=1,itochk
              if(ires_tot_tail(ito).gt.jstart) then
                itochk_tmp=itochk_tmp+1
                string66_tail_tmp(itochk_tmp)=string66_tail(ito)
                ires_tot_tail_tmp(itochk_tmp)=ires_tot_tail(ito)
                dtail_tmp(itochk_tmp)=dtail(ito)
                xtail_tmp(itochk_tmp)=xtail(ito)
                ytail_tmp(itochk_tmp)=ytail(ito)
                ztail_tmp(itochk_tmp)=ztail(ito)
              endif
              if(ires_tot_tail(ito).eq.jstart+1) then

C v4.5 since each residue might now have associated with it a bunch of
C HETATM entries we need to make sure that we're only looking at ATOMs

                if(string66_tail(ito)(1:4).ne.'ATOM') cycle

                if(string66_tail(ito)(13:16).eq.' N  ') then
                  rfx1=xtail(ito)
                  rfy1=ytail(ito)
                  rfz1=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' CA ') then
                  rfx2=xtail(ito)
                  rfy2=ytail(ito)
                  rfz2=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' C  ') then
                  rfx3=xtail(ito)
                  rfy3=ytail(ito)
                  rfz3=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' O  ') then
                  rfx4=xtail(ito)
                  rfy4=ytail(ito)
                  rfz4=ztail(ito)
                endif
              endif
            enddo
            do ito=1,itochk_tmp
              string66_tail(ito)=string66_tail_tmp(ito)
              ires_tot_tail(ito)=ires_tot_tail_tmp(ito)
              dtail(ito)=dtail_tmp(ito)
              xtail(ito)=xtail_tmp(ito)
              ytail(ito)=ytail_tmp(ito)
              ztail(ito)=ztail_tmp(ito)
            enddo
            itochk=itochk_tmp

C do a sanity check here

            if(itochk.eq.0) then
              if(jbacktrack.ne.iend_tot_nters(i)) then
                write(*,*)
                write(*,*)'something wrong here'
                write(*,*)'jbacktrack = ',jbacktrack
                write(*,*)'iend_tot_nters(i) = ',iend_tot_nters(i)
                write(*,*)
                stop
              endif
            endif

C do C-terminal tails here

          elseif(typ_tot_tails(itmp).eq.2) then

C v4.4 check if we've failed nbacktries times - if we have then we'll
C try to reset jbacktrack further back up the chain

            nbackmoves=nbackmoves+1
            if(nbackmoves.ge.nbacktries) then

C if we've repeatedly restarted it may be time to start from scratch

              nbackresets=nbackresets+1
              if(nbackresets.ge.nbacktries) then

C first check if we need to delete structure from the tail - note that
C we only do this if the logical "iworkoutwards" is invoked - if it is
C not then we do nothing and just restart the tail

                ngiveups=ngiveups+1
                if(ngiveups.ge.nbacktries.and.iworkoutwards) then
                  kwant=0
                  do k=ibeg_tot_cters(i),iend_tot_cters(i)
                    if(inot_made_yet(k).eq.0) then
                      kwant=k
                      exit
                    endif
                  enddo
                  if(kwant.ne.0) then
                    write(*,*)'giving up^2; remove struc frm res#',kwant
                    inot_made_yet(kwant)=1
                  endif
                  ngiveups=0
                  nfails_chn=0 ! reset this so we give it a fair shake
                  goto 5432
                else
                  write(*,*)'giving up; start from scratch # ',ngiveups
                  goto 5432
                endif

              endif

C v4.4 count back to find another backtrack position: jbacktrack - note
C that we only count flexible residues in the count...

              jbacktrack_old=jbacktrack
              imax_chn=jbacktrack_old  ! reset max res# found for C-ter builds
              l=0
              kwant=0
              do k=jbacktrack-1,ibeg_tot_cters(i),-1
                if(inot_made_yet(k).eq.1) then
                  l=l+1
                  if(l.eq.nbacktrack) then
                    kwant=k
                    exit
                  endif
                endif
              enddo
              if(kwant.ne.0) then
                jbacktrack=kwant
              else
                jbacktrack=ibeg_tot_cters(i)
              endif
              nbackmoves=0
              if(jbacktrack.ne.jbacktrack_old) then
              write(*,*)'NEW backtrack position for model# ',nmodels+1,
     &                  ' jbacktrack = ',jbacktrack
              endif
            endif

            jstart=jbacktrack

C now keep only those atoms for which j<jstart
C v4.5 this will need checking to see if it works if we have HETATM
C entries - I think that as written it won't - TEMP TEMP TEMP - so for
C now we should only try HETATM-association with nbacktrack=0
 
            itochk_tmp=0
            do ito=1,itochk
              if(ires_tot_tail(ito).lt.jstart) then
                itochk_tmp=itochk_tmp+1
                string66_tail_tmp(itochk_tmp)=string66_tail(ito)
                ires_tot_tail_tmp(itochk_tmp)=ires_tot_tail(ito)
                dtail_tmp(itochk_tmp)=dtail(ito)
                xtail_tmp(itochk_tmp)=xtail(ito)
                ytail_tmp(itochk_tmp)=ytail(ito)
                ztail_tmp(itochk_tmp)=ztail(ito)
              endif
              if(ires_tot_tail(ito).eq.jstart-1) then

C v4.5 since each residue might now have associated with it a bunch of
C HETATM entries we need to make sure that we're only looking at ATOMs

                if(string66_tail(ito)(1:4).ne.'ATOM') cycle

                if(string66_tail(ito)(13:16).eq.' N  ') then
                  rfx1=xtail(ito)
                  rfy1=ytail(ito)
                  rfz1=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' CA ') then
                  rfx2=xtail(ito)
                  rfy2=ytail(ito)
                  rfz2=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' C  ') then
                  rfx3=xtail(ito)
                  rfy3=ytail(ito)
                  rfz3=ztail(ito)
                elseif(string66_tail(ito)(13:16).eq.' O  ') then
                  rfx4=xtail(ito)
                  rfy4=ytail(ito)
                  rfz4=ztail(ito)
                endif
              endif
            enddo
            do ito=1,itochk_tmp
              string66_tail(ito)=string66_tail_tmp(ito)
              ires_tot_tail(ito)=ires_tot_tail_tmp(ito)
              dtail(ito)=dtail_tmp(ito)
              xtail(ito)=xtail_tmp(ito)
              ytail(ito)=ytail_tmp(ito)
              ztail(ito)=ztail_tmp(ito)
            enddo
            itochk=itochk_tmp

C do a sanity check here

            if(itochk.eq.0) then
              if(jbacktrack.ne.ibeg_tot_cters(i)) then
                write(*,*)
                write(*,*)'something wrong here'
                write(*,*)'jbacktrack = ',jbacktrack
                write(*,*)'ibeg_tot_cters(i) = ',ibeg_tot_cters(i)
                write(*,*)
                stop
              endif
            endif

C close the if over type of tail...

          endif

C if itrybacktrack is zero then we just restart the tail from scratch

        elseif(itrybacktrack.eq.0) then
          itochk=0
        endif

C now read dipeptide structures in loopy_inp.pdb if they exist - note
C that we only do this on the first time through for each tail...

        if(isentback.eq.0) then

C it's probably a good thing that this section is called only
C infrequently as it looks massively inefficient - let's try to keep
C track of how slow it might be by putting in simple writes - for now I
C don't think it's worth trying to speed this up but keep it in mind

          write(*,*)
          write(*,*)'starting to find structured dipeptides'

C first reset all temp4 entries to zero

          n_temp4=0
          temp4_N=0
          temp4_CA=0
          temp4_C=0
          temp4_O=0

C v2.1 set residue numbers correctly here - anything that depends on a
C correctly numbered residue from here on out will use "jx" not "j"

          if(typ_tot_tails(itmp).eq.1) then
            jx_beg=ibeg_tot_nters(i)-my_frst_res(my_frag)+
     &                                cur_offset(my_frag)+1
            jx_end=iend_tot_nters(i)-my_frst_res(my_frag)+
     &                                cur_offset(my_frag)+1
          elseif(typ_tot_tails(itmp).eq.2) then
            jx_beg=ibeg_tot_cters(i)-my_frst_res(my_frag)+
     &                                cur_offset(my_frag)+1
            jx_end=iend_tot_cters(i)-my_frst_res(my_frag)+
     &                                cur_offset(my_frag)+1
          endif

C now read the loopy_inp.pdb for backbone atoms of all residues
C these atoms will be used for superposition operations...

          open(unit=41,file='loopy_inp.pdb',status='unknown')

          ntmp=0

511       read(41,'(a66)',end=512)string66

c         ntmp=ntmp+1 ! for debug
c         call system('echo "HI THERE 391"')
c         write(*,*)'istatus 391 = ',istatus,' ntmp= ',ntmp
c         if(istatus.ne.0) stop ! TEMP TEMP TEMP

C v4.5 4 April 2021 - skip if it's a HETATM? - CRAZY I ONLY JUST FOUND

          if(string66(1:4).ne.'ATOM') goto 511

C following line is a BUGFIX added late in the day to v1.2

          if(string66(55:60).eq.'-99.99') goto 511 ! skip fake atoms!

C skip if it's the wrong chain ID

          if(string66(22:22).ne.cur_chainID(my_frag)) goto 511

C v2.5 - I think we should read all res's in this chain in case the tail
C extends later into this area - note that we only read this file
C *once* per tail so read all possible res's now...

          read(string66(23:26),'(i4)')ires

C n_temp4 is the # of atoms in this residue

          n_temp4(ires)=n_temp4(ires)+1
          string66_temp4(n_temp4(ires),ires)=string66 
          read(string66,'(30x,3f8.3)')x,y,z
          d_temp4(n_temp4(ires),ires)=0 ! v4.5 peptide default to dom=0
          x_temp4(n_temp4(ires),ires)=x
          y_temp4(n_temp4(ires),ires)=y
          z_temp4(n_temp4(ires),ires)=z
          if(string66(13:16).eq.' N  ') temp4_N(ires) =n_temp4(ires)
          if(string66(13:16).eq.' CA ') temp4_CA(ires)=n_temp4(ires)
          if(string66(13:16).eq.' C  ') temp4_C(ires) =n_temp4(ires)
          if(string66(13:16).eq.' O  ') temp4_O(ires) =n_temp4(ires)

          goto 511

512       close(41) ! v4.5 should be a close - trying debug sys calls
c512      rewind(12)

ccc       stop ! TEMP TEMP TEMP

C now loop over all residues in the tail and decide if we have a
C complete structured dipeptide for them or not...

          if(typ_tot_tails(itmp).eq.1) then

            do j=ibeg_tot_nters(i),iend_tot_nters(i)

              jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

C default to assuming that we have not yet made both residues (j,j+1) of
C this dipeptide - if inot_made_yet=1 later we will use our dipeptide
C libraries to sample structures, otherwise we use our temp4-like files

              inot_made_yet(j)=1

C if either of the two residues don't exist then we can cycle now

              if(temp4_CA(jx).eq.0.or.temp4_CA(jx+1).eq.0) cycle

C make sure that the CA atoms are not separated by crazy distances -
C this should only be necessary in cases where we're building a
C structure using premade domains (i.e. using "rebuild_tails" or
C "auto_rebuild")...

              iforget=-1
              x1=x_temp4(temp4_CA(jx),jx)
              y1=y_temp4(temp4_CA(jx),jx)
              z1=z_temp4(temp4_CA(jx),jx)
              x2=x_temp4(temp4_CA(jx+1),jx+1)
              y2=y_temp4(temp4_CA(jx+1),jx+1)
              z2=z_temp4(temp4_CA(jx+1),jx+1)

              dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
              if(dist.ge.5.0) then
                write(*,*)'crazy connection b/w N-ter struct res ',j,j+1
                write(*,*)'this will be replaced with dipeptide library'
                iforget=1
              else
                iforget=0
              endif

              if(temp4_N(jx).ne.0.and.temp4_N(jx+1).ne.0.and.
     &           iforget.eq.0) inot_made_yet(j)=0 ! note that we found this dipeptide

            enddo 

C v4.5 - at this point we can look through out HETATM list and add on
C any entries that should be associated with this residue - we look for
C the first structured residue and add the HETATM entries to it if they
C are within 500A of its CA atom

            mydom_last=-999

            write(*,*)'finding attached HETATMs for N-ter tail# ',i

            do j=iend_tot_nters(i),ibeg_tot_nters(i),-1

              if(inot_made_yet(j).eq.0) then

                mydom=i_in_dom_fasta_tot(j)
                write(*,*)'check me ',j,mydom

C if this structured residue is part of the same domain as the previous
C structured residue that we dealt with then skip it

                if(mydom.eq.mydom_last) cycle

C otherwise it's time to see if any HETATM entries are associated with
C this new domain...

                naddons=0
                jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

                x1=x_temp4(temp4_CA(jx),jx)
                y1=y_temp4(temp4_CA(jx),jx)
                z1=z_temp4(temp4_CA(jx),jx)

                do n=1,num_het_atms
                
                  if(ih(n).eq.0) cycle
 
                  x2=xh(n)
                  y2=yh(n)
                  z2=zh(n)
                  dist2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                  if(dist2.le.250000.0) then ! within 500A
                    naddons=naddons+1
                    n_temp4(jx)=n_temp4(jx)+1
                    d_temp4(n_temp4(jx),jx)       =mydom
                    x_temp4(n_temp4(jx),jx)       =xh(n)
                    y_temp4(n_temp4(jx),jx)       =yh(n)
                    z_temp4(n_temp4(jx),jx)       =zh(n)
                    string66_temp4(n_temp4(jx),jx)=sh(n) 

C record that we no longer want this HETATM from HETERO_original.pdb

                    ih(n)=0

                  endif
 
                enddo

C record this domain # in mydom_last so we don't do it again

                mydom_last=mydom

              endif
 
            enddo

            write(*,*)'found ',naddons,' HETATMs to add ',
     &                'to domain# ',mydom,' in N-ter tail # ',i

C do the same code for C-terminal tails but with jx-1, not jx+1

          elseif(typ_tot_tails(itmp).eq.2) then

            do j=ibeg_tot_cters(i),iend_tot_cters(i)

              jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

C default to assuming that we have not yet made both residues (j,j-1) of
C this dipeptide - if inot_made_yet=1 later we will use our dipeptide
C libraries to sample structures, otherwise we use our temp4-like files

              inot_made_yet(j)=1

C if either of the two residues don't exist then we can cycle now

              if(temp4_CA(jx).eq.0.or.temp4_CA(jx-1).eq.0) cycle

C make sure that the CA atoms are not separated by crazy distances -
C this should only be necessary in cases where we're building a
C structure using premade domains (i.e. using "rebuild_tails" or
C "auto_rebuild")...

              iforget=-1
              x1=x_temp4(temp4_CA(jx),jx)
              y1=y_temp4(temp4_CA(jx),jx)
              z1=z_temp4(temp4_CA(jx),jx)
              x2=x_temp4(temp4_CA(jx-1),jx-1)
              y2=y_temp4(temp4_CA(jx-1),jx-1)
              z2=z_temp4(temp4_CA(jx-1),jx-1)

              dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
              if(dist.ge.5.0) then
                write(*,*)'crazy connection b/w C-ter struct res ',j,j-1
                write(*,*)'this will be replaced with dipeptide library'
                iforget=1
              else
                iforget=0
              endif

              if(temp4_N(jx).ne.0.and.temp4_N(jx-1).ne.0.and.
     &           iforget.eq.0) inot_made_yet(j)=0 ! note that we found this dipeptide

            enddo 

C v4.5 - at this point we can look through out HETATM list and add on
C any entries that should be associated with this residue - we look for
C the first structured residue and add the HETATM entries to it if they
C are within 500A of its CA atom

            mydom_last=-999

            write(*,*)'finding attached HETATMs for C-ter tail# ',i

C v4.5REDUX - there was a bug here - till 18 March 2021 this said
C nters instead of cters...

            do j=ibeg_tot_cters(i),iend_tot_cters(i)

              if(inot_made_yet(j).eq.0) then

                mydom=i_in_dom_fasta_tot(j)

C if this structured residue is part of the same domain as the previous
C structured residue that we dealt with then skip it

                if(mydom.eq.mydom_last) cycle

C otherwise it's time to see if any HETATM entries are associated with
C this new domain...

                naddons=0
                jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

                x1=x_temp4(temp4_CA(jx),jx)
                y1=y_temp4(temp4_CA(jx),jx)
                z1=z_temp4(temp4_CA(jx),jx)

ccc             write(*,*)'JFC CHECK ',i,j,mydom,x1,y1,z1
   
                do n=1,num_het_atms
                 
                  if(ih(n).eq.0) cycle
 
                  x2=xh(n)
                  y2=yh(n)
                  z2=zh(n)
                  dist2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
                  if(dist2.le.250000.0) then ! within 500A
                    naddons=naddons+1
ccc     if(i.eq.2) then
ccc       write(*,8519)i,naddons,x1,y1,z1,n,ih(n),x2,y2,z2,sqrt(dist2)
8519      format('WTF ',2i8,3f12.3,2i8,4f12.3)
ccc     endif
                    n_temp4(jx)=n_temp4(jx)+1
                    d_temp4(n_temp4(jx),jx)       =mydom
                    x_temp4(n_temp4(jx),jx)       =xh(n)
                    y_temp4(n_temp4(jx),jx)       =yh(n)
                    z_temp4(n_temp4(jx),jx)       =zh(n)
                    string66_temp4(n_temp4(jx),jx)=sh(n) 

C record that we no longer want this HETATM from HETERO_original.pdb

                    ih(n)=0

                  endif
 
                enddo

C record this domain # in mydom_last so we don't do it again

                mydom_last=mydom

              endif
 
            enddo

            write(*,*)'found ',naddons,' HETATMs to add ',
     &                'to domain# ',mydom,' in C-ter tail # ',i

          endif

C need to remember to close(12)

ccc       close(12)

          write(*,*)'done with finding structured dipeptides'
          write(*,*)

        endif

C now start growing the N-terminal tail from the end to the beginning

        if(typ_tot_tails(itmp).eq.1) then

C v4.3

ccc       do j=iend_tot_nters(i),ibeg_tot_nters(i),-1
          do j=jstart,ibeg_tot_nters(i),-1

C v4.4 set appropriate value for frelax

            if(j.gt.jrelax) then
              clash_fac=frelax
              clash_fac2=clash_fac**2 ! needed cos we deal with dist2
ccc           write(*,*)'using relaxed clash criteria for res# ',j
            else
              clash_fac=1.0d0
              clash_fac2=1.0d0
            endif
c           if(nmodels.eq.0) then
c           write(*,*)'for N-ter res# ',j,' using clash_fac= ',clash_fac
c           endif

c           write(*,*)'trying to add res # ',j,' to Nter chain # ',i
            nfails_res=0 ! #failures for this res

            nt1=j-my_frst_res(ifrg_tot_nters(i))
            if(nt1.lt.imin_chn) then
              imin_chn=nt1

C v4.4 count up to find the new backtrack position: jbacktrack - note
C that we only count flexible residues in the count...

              l=0
              kwant=0
              do k=j+1,iend_tot_nters(i)
                if(inot_made_yet(k).eq.1) then
                  l=l+1
                  if(l.eq.nbacktrack) then
                    kwant=k
                    exit
                  endif
                endif
              enddo
              if(kwant.ne.0) then
                jbacktrack=kwant
              else
                jbacktrack=iend_tot_nters(i)
              endif
ccccc         nbackmoves=0
              write(*,*)'high-tide mark = ',j,' for model# ',nmodels+1,
     &                  ' jbacktrack = ',jbacktrack
            endif
 
C v2.1 set residue numbers correctly here - anything that depends on a
C correctly numbered residue from here on out will use "jx" not "j"

            jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

C v3.0 keep a chain-local residue number also for PSIPRED matching

            jloc=j-my_frst_res(my_frag)+1

C we now need the backbone coords of the residue after j - this will be
C used to place the new dipeptide... - first, decide which loopy pdb
C file to open (depends on if this is first res of the tail or not)

C v4.3 THIS WILL NEED CHANGING - IF THE FIRST BUT WE'VE DONE A READ-BACK
C THEN WE WON'T BE USING TEMP4 ANYMORE...
C v4.5REDUX TEMP TEMP TEMP

            if(j.eq.iend_tot_nters(i)) then
              rf(1,1)=x_temp4(temp4_N(jx+1),jx+1)
              rf(1,2)=y_temp4(temp4_N(jx+1),jx+1)
              rf(1,3)=z_temp4(temp4_N(jx+1),jx+1)
              rf(2,1)=x_temp4(temp4_CA(jx+1),jx+1)
              rf(2,2)=y_temp4(temp4_CA(jx+1),jx+1)
              rf(2,3)=z_temp4(temp4_CA(jx+1),jx+1)
              rf(3,1)=x_temp4(temp4_C(jx+1),jx+1)
              rf(3,2)=y_temp4(temp4_C(jx+1),jx+1)
              rf(3,3)=z_temp4(temp4_C(jx+1),jx+1)
              rf(4,1)=x_temp4(temp4_O(jx+1),jx+1)
              rf(4,2)=y_temp4(temp4_O(jx+1),jx+1)
              rf(4,3)=z_temp4(temp4_O(jx+1),jx+1)
            elseif(j.eq.jstart) then
              rf(1,1)=rfx1
              rf(1,2)=rfy1
              rf(1,3)=rfz1
              rf(2,1)=rfx2
              rf(2,2)=rfy2
              rf(2,3)=rfz2
              rf(3,1)=rfx3
              rf(3,2)=rfy3
              rf(3,3)=rfz3
              rf(4,1)=rfx4
              rf(4,2)=rfy4
              rf(4,3)=rfz4
            else
              rf(1,1)=rf_prev(1,1)
              rf(1,2)=rf_prev(1,2)
              rf(1,3)=rf_prev(1,3)
              rf(2,1)=rf_prev(2,1)
              rf(2,2)=rf_prev(2,2)
              rf(2,3)=rf_prev(2,3)
              rf(3,1)=rf_prev(3,1)
              rf(3,2)=rf_prev(3,2)
              rf(3,3)=rf_prev(3,3)
              rf(4,1)=rf_prev(4,1)
              rf(4,2)=rf_prev(4,2)
              rf(4,3)=rf_prev(4,3)
            endif

C get ready to select a dipeptide conformation - here we're just
C preparing to figure out the filename...

            dipep_cur(1:1)=seq_tot_res(j:j)
            dipep_cur(2:2)=seq_tot_res(j+1:j+1)
            do k=1,20
              if(dipep_cur(1:1).eq.dipep(k:k)) i1=k
              if(dipep_cur(2:2).eq.dipep(k:k)) i2=k
            enddo

C we need a fast way of selecting conformations without repeatedly
C selecting the same residue - to do this, generate a list of random
C numbers, find the min value and select that conf, then set its value
C to 2 so that it isn't selected again...

            do k=1,num_peps(i1,i2)
              call random_number(u1)
              rand_res(k)=u1
            enddo

C come back here if we failed to place this conformation - note that we
C had to pull this statement out of the inot_made_yet if statment...

808         r=minval(rand_res(1:num_peps(i1,i2)))
 
C now we need to determine which dipeptide pdb file to read - we could
C either be reading from a library or we could be reading from pre-built
C do the following if reading from the dipeptide library...

            if(inot_made_yet(j).eq.1) then

C make sure that we set i_am_original for this residue to zero...

              i_am_original(j)=0

c             if(nfails_res.eq.0) 
c    &        write(*,*)'reading library for res# ',j,j+1,dipep_cur

C v3.0 randomly select whether to put in a helical residue here
C basic form is correct here but we need to be sure my_PSIPRED_vl2 
C is accessed correctly - jkl should be sequence# corresponding to this
C fragment, 1 should always be correct (we only allow one pdb to act as
C the template), and jloc should be the chain-local res# - so it all
C looks correct at this stage but obviously needs testing...

C v3.4 instead we will automatically select helix if the confidence
C value is above the cutoff specified in the input file:
C   helix_cut

C v3.5 implement same idea for sheets...

C note that we don't allow for dipeptides to have HE or EH conformations
C (i.e. helix followed immediately by a sheet) - so, as written, the if
C statement below will give preference to helices: a HE segment will be
C treated as a HC segment, and a EH segment will become a CH segment

c             call random_number(v1)
              jkl=my_frags_seq_num(my_frag)
c             write(*,*)'check PSIPRED ',jloc,
c    &                      my_PSIPRED_vl2(jkl,1,jloc)

Cv3.0         if(v1.lt.my_PSIPRED_vl2(jkl,1,jloc)) then
Cv3.0           k=1
Cv3.0         else
Cv3.0           k=minloc(rand_res(1:num_peps(i1,i2)),1)
Cv3.0           rand_res(k)=2.0 ! make sure k is not selected again
Cv3.0         endif

              if(my_PSIPRED_vl2(jkl,1,jloc+1).lt.helix_cut.and.
     &           my_PSIPRED_vl2(jkl,1,jloc  ).lt.helix_cut.and.
     &           my_PSIPRED_vl3(jkl,1,jloc+1).lt.sheet_cut.and.
     &           my_PSIPRED_vl3(jkl,1,jloc  ).lt.sheet_cut) then
                k=minloc(rand_res(1:dipep_both_coils_end(i1,i2)),1)
              elseif(my_PSIPRED_vl2(jkl,1,jloc+1).ge.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).lt.helix_cut) then
                k=minloc(rand_res(dipep_helix_Nter_beg(i1,i2):
     &                            dipep_helix_Nter_end(i1,i2)),1)
                k=k+dipep_helix_Nter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl2(jkl,1,jloc+1).ge.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).ge.helix_cut) then
                k=minloc(rand_res(dipep_both_helix_beg(i1,i2):
     &                            dipep_both_helix_end(i1,i2)),1)
                k=k+dipep_both_helix_beg(i1,i2)-1
              elseif(my_PSIPRED_vl2(jkl,1,jloc+1).lt.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).ge.helix_cut) then
                k=minloc(rand_res(dipep_helix_Cter_beg(i1,i2):
     &                            dipep_helix_Cter_end(i1,i2)),1)
                k=k+dipep_helix_Cter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc+1).ge.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).lt.sheet_cut) then
                k=minloc(rand_res(dipep_sheet_Nter_beg(i1,i2):
     &                            dipep_sheet_Nter_end(i1,i2)),1)
                k=k+dipep_sheet_Nter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc+1).ge.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).ge.sheet_cut) then
                k=minloc(rand_res(dipep_both_sheet_beg(i1,i2):
     &                            dipep_both_sheet_end(i1,i2)),1)
                k=k+dipep_both_sheet_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc+1).lt.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).ge.sheet_cut) then
                k=minloc(rand_res(dipep_sheet_Cter_beg(i1,i2):
     &                            dipep_sheet_Cter_end(i1,i2)),1)
                k=k+dipep_sheet_Cter_beg(i1,i2)-1
              endif

c         write(*,47816)jloc,jloc+1,
c    &                        my_PSIPRED_vl2(jkl,1,jloc),
c    &                        my_PSIPRED_vl2(jkl,1,jloc+1),
c    &                        dipep_both_coils_beg(i1,i2),
c    &                        dipep_both_coils_end(i1,i2),
c    &                        dipep_helix_Nter_beg(i1,i2),
c    &                        dipep_helix_Nter_end(i1,i2),
c    &                        dipep_both_helix_beg(i1,i2),
c    &                        dipep_both_helix_end(i1,i2),
c    &                        dipep_helix_Cter_beg(i1,i2),
c    &                        dipep_helix_Cter_end(i1,i2),k,
c    &                        num_peps(i1,i2)
47816       format('checkage ',2i6,2f8.3,11i6)

              rand_res(k)=2.0 ! make sure k is not selected again

              xm=0.0 ! remember to zero out all entries
              ym=0.0 ! remember to zero out all entries
              zm=0.0 ! remember to zero out all entries
              dm=0   ! v4.5 remember to zero out all entries
              xm(1)=x_dipep(dipep_N2(i1,i2),k,i1,i2)
              ym(1)=y_dipep(dipep_N2(i1,i2),k,i1,i2)
              zm(1)=z_dipep(dipep_N2(i1,i2),k,i1,i2)
              xm(2)=x_dipep(dipep_CA2(i1,i2),k,i1,i2)
              ym(2)=y_dipep(dipep_CA2(i1,i2),k,i1,i2)
              zm(2)=z_dipep(dipep_CA2(i1,i2),k,i1,i2)
              xm(3)=x_dipep(dipep_C2(i1,i2),k,i1,i2)
              ym(3)=y_dipep(dipep_C2(i1,i2),k,i1,i2)
              zm(3)=z_dipep(dipep_C2(i1,i2),k,i1,i2)
              xm(4)=x_dipep(dipep_O2(i1,i2),k,i1,i2)
              ym(4)=y_dipep(dipep_O2(i1,i2),k,i1,i2)
              zm(4)=z_dipep(dipep_O2(i1,i2),k,i1,i2)

c      write(*,*)'BEF1',xm(1),ym(1),zm(1),rf(1,1),rf(1,2),rf(1,3)
c      write(*,*)'BEF2',xm(2),ym(2),zm(2),rf(2,1),rf(2,2),rf(2,3)
c      write(*,*)'BEF3',xm(3),ym(3),zm(3),rf(3,1),rf(3,2),rf(3,3)
c      write(*,*)'BEF4',xm(4),ym(4),zm(4),rf(4,1),rf(4,2),rf(4,3)

              do l=1,dipep_natom1(i1,i2)
                xm(4+l)=x_dipep(l,k,i1,i2)
                ym(4+l)=y_dipep(l,k,i1,i2)
                zm(4+l)=z_dipep(l,k,i1,i2)
              enddo  
              nm=4+dipep_natom1(i1,i2)
              itoadd=dipep_natom1(i1,i2) ! store how many atoms we will add
              itoskp=0                   ! store how many atoms to skip

            elseif(inot_made_yet(j).eq.0) then

C do the following if reading the dipeptide from pre-built structure
C note that we have a simple way to switch back to use the dipeptide
C library if we fail excessively with the pre-built structure: we just
C reset inot_made_yet to 1 for this particular dipeptide...

c             if(nfails_res.eq.0) 
c    &        write(*,*)'reading strctre for res# ',j,j+1,dipep_cur

              xm=0.0 ! remember to zero out all entries
              ym=0.0 ! remember to zero out all entries
              zm=0.0 ! remember to zero out all entries
              dm=0   ! v4.5 remember to zero out all entries
              xm(1)=x_temp4(temp4_N(jx+1),jx+1)
              ym(1)=y_temp4(temp4_N(jx+1),jx+1)
              zm(1)=z_temp4(temp4_N(jx+1),jx+1)
              xm(2)=x_temp4(temp4_CA(jx+1),jx+1)
              ym(2)=y_temp4(temp4_CA(jx+1),jx+1)
              zm(2)=z_temp4(temp4_CA(jx+1),jx+1)
              xm(3)=x_temp4(temp4_C(jx+1),jx+1)
              ym(3)=y_temp4(temp4_C(jx+1),jx+1)
              zm(3)=z_temp4(temp4_C(jx+1),jx+1)
              xm(4)=x_temp4(temp4_O(jx+1),jx+1)
              ym(4)=y_temp4(temp4_O(jx+1),jx+1)
              zm(4)=z_temp4(temp4_O(jx+1),jx+1)
              do l=1,n_temp4(jx)
                xm(4+l)=x_temp4(l,jx)
                ym(4+l)=y_temp4(l,jx)
                zm(4+l)=z_temp4(l,jx)
                dm(4+l)=d_temp4(l,jx) ! v4.5
              enddo  
              nm=4+n_temp4(jx)
              itoadd=n_temp4(jx) ! store how many atoms we'll add
              itoskp=0           ! store how many atoms to skip

            endif

C at this point we should have our files 'target' and 'temp4' and be
C ready to do a superposition...

            call pdbsup(nm,xm,ym,zm,rf,rms)
c           write(*,*)'rms = ',rms,' nm = ',nm
c           write(91,'(f10.5)')rms

c      write(*,*)'AFT1',xm(1),ym(1),zm(1),rf(1,1),rf(1,2),rf(1,3)
c      write(*,*)'AFT2',xm(2),ym(2),zm(2),rf(2,1),rf(2,2),rf(2,3)
c      write(*,*)'AFT3',xm(3),ym(3),zm(3),rf(3,1),rf(3,2),rf(3,3)
c      write(*,*)'AFT4',xm(4),ym(4),zm(4),rf(4,1),rf(4,2),rf(4,3)

C now decide whether to accept this on the basis of its rms...

            if(rms.gt.rms_crit) then

              nfails_res=nfails_res+1
ccc           write(*,*)'HELP1N!! ',nfails_res

C note that the last criterion ensures that we immediately quit if a
C structured peptide was used to add the new residue

              if(nfails_res.ge.num_peps(i1,i2).or.
     &           nfails_res.ge.nfails_res_max.or.
     &           inot_made_yet(j).eq.0) then

                nfails_chn=nfails_chn+1
 
C 2023_v1.1 - only extend if nfails_chn_max>0
 
                if(nfails_chn.ge.nfails_chn_max.and.
     &             nfails_chn_max.gt.0) then
                  write(*,*)'rms   extending Nter chain # ',i,rms
                  isentback=1
                  goto 404 ! extend the chain
                else
c                 write(*,*)'rms   restartng Nter chain # ',i,rms

C v4.3 make sure we try to backtrack

                  if(nbacktrack.gt.0) itrybacktrack=1

                  isentback=1 ! BUGFIX ADDED LATE TO SPEED UP
                  goto 505 ! restart the chain
                endif
  
              else
                goto 808 ! pick another conformation
              endif

            endif

C if we got here then we may be about to include the new amino acid...
C first, get the atoms that we want from the superimpose output

            do m=1,itoadd
              xm(m)=xm(m+4)
              ym(m)=ym(m+4)
              zm(m)=zm(m+4)
              dm(m)=dm(m+4) ! v4.5
              iresm(m)=jx ! v2.1 store res# for debugging
            enddo

C if membrane_dive then get the z coords of the N atoms

            if(membrane_dive) then
              z1=zm(1+4) ! store z coord of N atom in new res
     &                   ! this is always first atom
              z2=rf(1,3) ! store z coord of 2nd N atom - 1st fixd atom
ccc           z2=zm(n_temp4(jx)+1+4) ! store z coord of 2nd N atom
            endif

C if membrane_dive and this is the first residue we're adding then we
C can use the previous residue to tell us which zone we're in - we've
C got four possibilities:

C (1) i_go_up=1 : tail starts below membrane (must stay that side)
C (2) i_go_up=2 : tail starts within membrane (and must move in neg)
C (3) i_go_up=3 : tail starts within membrane (and must move in pos)
C (4) i_go_up=4 : tail starts above membrane (must stay that side)

            if(membrane_dive.and.j.eq.iend_tot_nters(i)) then
              if(zmem_neg2.eq.0.000) then
                if(z2.le.zmem_neg1) then
                  i_go_up=1
                elseif(z2.le.0.0.and.z2.gt.zmem_neg1) then
                  i_go_up=2
                elseif(z2.le.zmem_pos1.and.z2.gt.0.0) then
                  i_go_up=3
                else
                  i_go_up=4
                endif
              else
                if(z2.le.zmem_neg1) then
                  i_go_up=1
                elseif(z2.le.zmem_ave1.and.z2.gt.zmem_neg1) then
                  i_go_up=2
                elseif(z2.le.zmem_pos1.and.z2.gt.zmem_ave1) then
                  i_go_up=3
                elseif(z2.le.zmem_neg2.and.z2.gt.zmem_pos1) then
                  i_go_up=4 
                elseif(z2.le.zmem_ave2.and.z2.gt.zmem_neg2) then
                  i_go_up=5
                elseif(z2.le.zmem_pos2.and.z2.gt.zmem_ave2) then
                  i_go_up=6
                else
                  i_go_up=7
                endif
              endif
            endif
        
C now, if membrane_dive, we can ask if we're going in right direction
C note that imem_quit - which is used for quitting this tail conf - is
C only set non-zero  if membrane_dive is true... remember that z1 is the
C new residue and z2 is the old residue...

            imem_quit=0 
            if(membrane_dive) then
              if(zmem_neg2.eq.0.000) then
                if(i_go_up.eq.1) then
                  if(z1.gt.zmem_neg1) imem_quit=1
                elseif(i_go_up.eq.2) then
                  if(z2.ge.zmem_neg1.and.z1.gt.z2) imem_quit=1
                elseif(i_go_up.eq.3) then
                  if(z2.le.zmem_pos1.and.z1.lt.z2) imem_quit=1
                elseif(i_go_up.eq.4) then
                  if(z1.lt.zmem_pos1) imem_quit=1
                endif
              else
                if(i_go_up.eq.1) then
                  if(z1.gt.zmem_neg1) imem_quit=1
                elseif(i_go_up.eq.2) then
                  if(z2.ge.zmem_neg1.and.z1.gt.z2) imem_quit=1
                elseif(i_go_up.eq.3) then
                  if(z2.le.zmem_pos1.and.z1.lt.z2) imem_quit=1
                elseif(i_go_up.eq.4) then
                  if(z1.lt.zmem_pos1) imem_quit=1
                  if(z1.gt.zmem_neg2) imem_quit=1
                elseif(i_go_up.eq.5) then
                  if(z2.ge.zmem_neg2.and.z1.gt.z2) imem_quit=1
                elseif(i_go_up.eq.6) then
                  if(z2.le.zmem_pos2.and.z1.lt.z2) imem_quit=1
                elseif(i_go_up.eq.7) then
                  if(z1.lt.zmem_pos2) imem_quit=1
                endif
              endif
c             write(*,*)
c             write(*,*)'check z1 z2 Nter settings ',z1,z2,i_go_up
c             write(*,*)
            endif
 
C now, before finally adding this residue we need to do a clash check -
C we'll use the same code to implement the quit due to membrane_dive
C rather than adding in all the same code above

C TEMP TEMP TEMP we need to remember that we are currently only checking
C against fixed residues - this will be everything except HETATM entries
C and except for other residues in the tail....

            iclash=0
            do m1=1,itoadd
              xt=xm(m1)
              yt=ym(m1)
              zt=zm(m1)

C skip clash-checking if this atom is not on the grid - note that this
C used to be a fatal error when I made a smaller clash grid that only
C encompassed the atoms near to the beginning of the tail...

              if(xt.lt.xmin.or.xt.gt.xmax.or.
     &           yt.lt.ymin.or.yt.gt.ymax.or.
     &           zt.lt.zmin.or.zt.gt.zmax) cycle

              it1=int((xt-xmin)*xinv)+1
              jt1=int((yt-ymin)*yinv)+1
              kt1=int((zt-zmin)*zinv)+1
              do m2=1,temp_grid(it1,jt1,kt1)

C skip the grid atom if it's the jx+1'th residue of same chain

                if(cchn_grid(m2,it1,jt1,kt1).eq.
     &             cur_chainID(my_frag).and.
     &             rchn_grid(m2,it1,jt1,kt1).eq.jx+1) then
c                 write(*,*)'skip clash with ',rchn_grid(m2,it1,jt1,kt1),
c    &                                         cchn_grid(m2,it1,jt1,kt1)
                  cycle
                endif

C otherwise, measure the distance and look for a clash

                xu=xatm_grid(m2,it1,jt1,kt1)
                yu=yatm_grid(m2,it1,jt1,kt1)
                zu=zatm_grid(m2,it1,jt1,kt1)
                dist2=((xt-xu)**2+(yt-yu)**2+(zt-zu)**2)

C v4.2 do this if it's a potential ATOM-ATOM clash 

                if(ihet_grid(m2,it1,jt1,kt1).eq.0) then
                  if(dist2.lt.clash_aa2*clash_fac2) then
                    iclash=1
                    jclash=rchn_grid(m2,it1,jt1,kt1)
c                   write(*,*)'intermolec clash ',iresm(m1),
c    &                rchn_grid(m2,it1,jt1,kt1),sqrt(dist2)
                    goto 9051
                  endif

C v4.2 do this if it's a potential ATOM-HETATM clash 

                elseif(ihet_grid(m2,it1,jt1,kt1).eq.1) then
                  if(dist2.lt.clash_ha2*clash_fac2) then
                    iclash=1
                    jclash=rchn_grid(m2,it1,jt1,kt1)
c                   write(*,*)'intermolec clash ',iresm(m1),
c    &                rchn_grid(m2,it1,jt1,kt1),sqrt(dist2)
                    goto 9051
                  endif
                endif

              enddo
            enddo

C check for intra-tail clashes - note that itochk contains residues j+1
C upwards so we need to explicitly check to ignore j:j+1 clashes

            do m1=1,itoadd
              xt=xm(m1)
              yt=ym(m1)
              zt=zm(m1)
              do m2=1,itochk
                xu=xtail(m2)
                yu=ytail(m2)
                zu=ztail(m2)
                dist2=((xt-xu)**2+(yt-yu)**2+(zt-zu)**2)
C v4.2          if(dist2.lt.clash2) then
                if(dist2.lt.clash_aa2*clash_fac2) then
                  iclash=1

C v2.7 here we ignore clashes between adjacent residues

                  jclash=ires_tot_tail(m2)
                  if(abs(jclash-j).eq.1) then
                    iclash=0
                    cycle
                  endif

C we have to be careful not to quit if the apparent "clash" is between
C two residues that are part of the same structured domain...
C here, remember that "j" and "jclash" are both in total-residue count
C numbering so we can plug them directly in and find their domains -
C note also that they must, by definition, be part of the same molecule

C v4.4 - we also need to check that they still have the inot_made_yet
C status set to 0 - these could have been updated to 1 if we failed to
C build the tail many times... - add this extra condition to the if
C statement so we only come here if all these things are true:

c                 if(i_in_dom_fasta_tot(j).eq.
c    &               i_in_dom_fasta_tot(jclash).and.
c    &               i_in_dom_fasta_tot(j).gt.0) then

                  if(i_in_dom_fasta_tot(j).eq.
     &               i_in_dom_fasta_tot(jclash).and.
     &               i_in_dom_fasta_tot(j).gt.0.and.
     &               inot_made_yet(j).eq.0.and.
     &               inot_made_yet(jclash).eq.0) then

                    iclash=0 ! set back to zero to ignore the clash
                    cycle
c                   write(*,*)'ignoring intra-tail clash ',
c    &                         iresm(m1),iresn(m2),sqrt(dist2),
c    &                         i_in_dom_fasta_tot(j)


C v4.5 we must do something similar for ATOM-HETATM entries where the
C HETATM is associated with the tail... TEMP TEMP TEMP

                  elseif(i_in_dom_fasta_tot(j).eq.dtail(m2).and.
     &                   inot_made_yet(j).eq.0) then
                    iclash=0
                    cycle

                  else
c                   write(*,*)'deadly intra-tail clash ',
c    &                         iresm(m1),iresn(m2),sqrt(dist2)
c                   write(*,*)'xm,ym,zm = ',xt,yt,zt
c                   write(*,*)'xn,yn,zn = ',xu,yu,zu
                    goto 9051
                  endif
                endif
              enddo
            enddo

9051        continue ! come here if we had a clash...

            if(iclash.eq.1.or.imem_quit.eq.1) then

C v4.7 6/10/2022 if imem_quit=1 let's restart the chain no questions
C asked - if we don't do this then we might keep snaking along the
C membrane leading to a highly biased distribution...
C UPDATE - the code seems to hang if we do this :(

ccc           if(imem_quit.eq.1) then
ccc             goto 5432   ! restart the chain
ccc           endif

              nfails_res=nfails_res+1
ccc           write(*,*)'HELP2N!! ',nfails_res

              if(nfails_res.ge.num_peps(i1,i2).or.
     &           nfails_res.ge.nfails_res_max.or.
     &           inot_made_yet(j).eq.0) then

                nfails_chn=nfails_chn+1
 
C 2023_v1.1 - only extend if nfails_chn_max>0
 
                if(nfails_chn.ge.nfails_chn_max.and.
     &             nfails_chn_max.gt.0) then
                  write(*,*)'clash extending Nter chain # ',i,j,jclash,
     &                     sqrt(dist2)
                  isentback=1
                  goto 404 ! extend the chain
                else
c                 write(*,*)'clash restartng Nter chain # ',i,j,jclash,
c    &                       sqrt(dist2),nfails_chn

C v4.3 make sure we try to backtrack

                  if(nbacktrack.gt.0) itrybacktrack=1

                  isentback=1 ! BUGFIX ADDED LATE TO SPEED UP
                  goto 505 ! restart the chain
                endif

              else
                goto 808 ! pick another conformation
              endif

            endif ! done with the "if clashing..."
  
C if we got here then we really are adding this residue...
C now add the atoms of the new residue to string66_tail...

            if(inot_made_yet(j).eq.1) then
              do l=1,itoadd
                if(string66_dipep(l,k,i1,i2)(13:16).eq.' N  ') then
                  rf_prev(1,1)=xm(l)
                  rf_prev(1,2)=ym(l)
                  rf_prev(1,3)=zm(l)
                elseif(string66_dipep(l,k,i1,i2)(13:16).eq.' CA ') then
                  rf_prev(2,1)=xm(l)
                  rf_prev(2,2)=ym(l)
                  rf_prev(2,3)=zm(l)
                elseif(string66_dipep(l,k,i1,i2)(13:16).eq.' C  ') then
                  rf_prev(3,1)=xm(l)
                  rf_prev(3,2)=ym(l)
                  rf_prev(3,3)=zm(l)
                elseif(string66_dipep(l,k,i1,i2)(13:16).eq.' O  ') then
                  rf_prev(4,1)=xm(l)
                  rf_prev(4,2)=ym(l)
                  rf_prev(4,3)=zm(l)
                endif

C v2.7 store information in internal strings to avoid constant writing

                write(string66_tail_tmp(l),81)
     &                      string66_dipep(l,k,i1,i2)(1:21),
     &                      cur_chainID(my_frag),
     &                      jx,xm(l),ym(l),zm(l),
     &                      string66_dipep(l,k,i1,i2)(55:66)       ! v2.7
                ires_tot_tail_tmp(l)=j
                dtail_tmp(l)=dm(l)
                xtail_tmp(l)=xm(l)
                ytail_tmp(l)=ym(l)
                ztail_tmp(l)=zm(l)
81              format(a21,a,i4,4x,3f8.3,a12)

              enddo
            else
              do l=1,itoadd

C v4.5 since each residue might now have associated with it a bunch of
C HETATM entries we need to make sure that we're only looking at ATOMs

                if(string66_temp4(l,jx)(1:4).eq.'ATOM') then

                  if(string66_temp4(l,jx)(13:16).eq.' N  ') then
                    rf_prev(1,1)=xm(l)
                    rf_prev(1,2)=ym(l)
                    rf_prev(1,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' CA ') then
                    rf_prev(2,1)=xm(l)
                    rf_prev(2,2)=ym(l)
                    rf_prev(2,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' C  ') then
                    rf_prev(3,1)=xm(l)
                    rf_prev(3,2)=ym(l)
                    rf_prev(3,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' O  ') then
                    rf_prev(4,1)=xm(l)
                    rf_prev(4,2)=ym(l)
                    rf_prev(4,3)=zm(l)
                  endif

                endif

C v2.7 store information in internal strings to avoid constant writing
C v4.5 make sure we don't change the residue number of HETATMs

                if(string66_temp4(l,jx)(1:4).eq.'ATOM') then

                  write(string66_tail_tmp(l),81)
     &                      string66_temp4(l,jx)(1:21),     
     &                      cur_chainID(my_frag),
     &                      jx,xm(l),ym(l),zm(l),
     &                      string66_temp4(l,jx)(55:66)            ! v2.7

                elseif(string66_temp4(l,jx)(1:4).ne.'ATOM') then

                  write(string66_tail_tmp(l),94)
     &                      string66_temp4(l,jx)(1:21),     
ccc  &                      cur_chainID(my_frag), ! v4.6 make sure
     &                      fin_chainID(my_frag), ! chainID is final
     &                      string66_temp4(l,jx)(23:26),     
     &                      xm(l),ym(l),zm(l),
     &                      string66_temp4(l,jx)(55:66)            ! v2.7

94                format(a21,a,a4,4x,3f8.3,a12)

                endif

                ires_tot_tail_tmp(l)=j
                dtail_tmp(l)=dm(l)
                xtail_tmp(l)=xm(l)
                ytail_tmp(l)=ym(l)
                ztail_tmp(l)=zm(l)

              enddo
            endif

C v2.7 here we should put the new atoms on the grid for clash checking




C now, add on the other residues from current_tail.pdb...

            do ito=1,itochk
              string66_tail_tmp(itoadd+ito)=string66_tail(ito)
              ires_tot_tail_tmp(itoadd+ito)=ires_tot_tail(ito)
              dtail_tmp(itoadd+ito)=dtail(ito)
              xtail_tmp(itoadd+ito)=xtail(ito)
              ytail_tmp(itoadd+ito)=ytail(ito)
              ztail_tmp(itoadd+ito)=ztail(ito)
            enddo

C update itochk to include the new residue

            itochk=itochk+itoadd

C finally, copy back to string66_tail...

            do ito=1,itochk
              string66_tail(ito)=string66_tail_tmp(ito)
              ires_tot_tail(ito)=ires_tot_tail_tmp(ito)
              dtail(ito)=dtail_tmp(ito)
              xtail(ito)=xtail_tmp(ito)
              ytail(ito)=ytail_tmp(ito)
              ztail(ito)=ztail_tmp(ito)
            enddo

C if this is the final residue then we write out to current_tail.pdb

            if(j.eq.ibeg_tot_nters(i)) then

              open(unit=12,file='current_tail.pdb',status='unknown')
              do l=1,itochk
                write(12,'(a66)')string66_tail(l)
              enddo
              close(12)

            endif

          enddo

        elseif(typ_tot_tails(itmp).eq.2) then

C v4.5REDUX

ccc       do j=ibeg_tot_cters(i),iend_tot_cters(i)
          do j=jstart,iend_tot_cters(i)

C v4.5REDUX set appropriate value for frelax

            if(j.lt.jrelax) then
              clash_fac=frelax
              clash_fac2=clash_fac**2 ! needed cos we deal with dist2
ccc           write(*,*)'using relaxed clash criteria for res# ',j
            else
              clash_fac=1.0d0
              clash_fac2=1.0d0
            endif
c           if(nmodels.eq.0) then
c           write(*,*)'for C-ter res# ',j,' using clash_fac= ',clash_fac
c           endif

c           write(*,*)'trying to add res # ',j,' to Cter chain # ',i
            nfails_res=0 ! #failures for this res

C v4.5REDUX - try to make this part of code echo the Nter stuff from 755
C lines or so earlier...

c           nt1=j-my_frst_res(ifrg_tot_cters(i))
c           if(nt1.gt.imax_chn) then
c             write(*,*)'high-tide mark = ',j,' for model# ',nmodels+1
c             imax_chn=nt1
c           endif
 
            nt1=j-my_frst_res(ifrg_tot_cters(i))
            if(nt1.gt.imax_chn) then
              imax_chn=nt1

C v4.4 count down to find the new backtrack position: jbacktrack - note
C that we only count flexible residues in the count...

              l=0
              kwant=0
              do k=j-1,ibeg_tot_cters(i),-1
                if(inot_made_yet(k).eq.1) then
                  l=l+1
                  if(l.eq.nbacktrack) then
                    kwant=k
                    exit
                  endif
                endif
              enddo
              if(kwant.ne.0) then
                jbacktrack=kwant
              else
                jbacktrack=ibeg_tot_cters(i)
              endif
ccccc         nbackmoves=0
              write(*,*)'high-tide mark = ',j,' for model# ',nmodels+1,
     &                  ' jbacktrack = ',jbacktrack
            endif
 
C v2.1 set residue numbers correctly here - anything that depends on a
C correctly numbered residue from here on out will use "jx" not "j"

            jx=j-my_frst_res(my_frag)+cur_offset(my_frag)+1

C v3.0 keep a chain-local residue number also for PSIPRED matching

            jloc=j-my_frst_res(my_frag)+1

C we now need the backbone coords of the residue before j - this will be
C used to place the new dipeptide... - first, decide which loopy pdb
C file to open (depends on if this is first res of the tail or not)

C v4.5REDUX again try to ape the Nter code here
C v4.3 THIS WILL NEED CHANGING - IF THE FIRST BUT WE'VE DONE A READ-BACK
C THEN WE WON'T BE USING TEMP4 ANYMORE...
C v4.5REDUX TEMP TEMP TEMP

ccc         if(j.eq.ibeg_tot_cters_orig) then
            if(j.eq.ibeg_tot_cters(i)) then
              rf(1,1)=x_temp4(temp4_N(jx-1),jx-1)
              rf(1,2)=y_temp4(temp4_N(jx-1),jx-1)
              rf(1,3)=z_temp4(temp4_N(jx-1),jx-1)
              rf(2,1)=x_temp4(temp4_CA(jx-1),jx-1)
              rf(2,2)=y_temp4(temp4_CA(jx-1),jx-1)
              rf(2,3)=z_temp4(temp4_CA(jx-1),jx-1)
              rf(3,1)=x_temp4(temp4_C(jx-1),jx-1)
              rf(3,2)=y_temp4(temp4_C(jx-1),jx-1)
              rf(3,3)=z_temp4(temp4_C(jx-1),jx-1)
              rf(4,1)=x_temp4(temp4_O(jx-1),jx-1)
              rf(4,2)=y_temp4(temp4_O(jx-1),jx-1)
              rf(4,3)=z_temp4(temp4_O(jx-1),jx-1)
ccc         elseif(j.eq.ibeg_tot_cters(i)) then
            elseif(j.eq.jstart) then
              rf(1,1)=rfx1
              rf(1,2)=rfy1
              rf(1,3)=rfz1
              rf(2,1)=rfx2
              rf(2,2)=rfy2
              rf(2,3)=rfz2
              rf(3,1)=rfx3
              rf(3,2)=rfy3
              rf(3,3)=rfz3
              rf(4,1)=rfx4
              rf(4,2)=rfy4
              rf(4,3)=rfz4
            else
              rf(1,1)=rf_prev(1,1)
              rf(1,2)=rf_prev(1,2)
              rf(1,3)=rf_prev(1,3)
              rf(2,1)=rf_prev(2,1)
              rf(2,2)=rf_prev(2,2)
              rf(2,3)=rf_prev(2,3)
              rf(3,1)=rf_prev(3,1)
              rf(3,2)=rf_prev(3,2)
              rf(3,3)=rf_prev(3,3)
              rf(4,1)=rf_prev(4,1)
              rf(4,2)=rf_prev(4,2)
              rf(4,3)=rf_prev(4,3)
            endif

C get ready to select a dipeptide conformation - here we're just
C preparing to figure out the filename...

            dipep_cur(1:1)=seq_tot_res(j-1:j-1)
            dipep_cur(2:2)=seq_tot_res(j:j)
            do k=1,20
              if(dipep_cur(1:1).eq.dipep(k:k)) i1=k
              if(dipep_cur(2:2).eq.dipep(k:k)) i2=k
            enddo

C we need a fast way of selecting conformations without repeatedly
C selecting the same residue - to do this, generate a list of random
C numbers, find the min value and select that conf, then set its value
C to 2 so that it isn't selected again...

            do k=1,num_peps(i1,i2)
              call random_number(u1)
              rand_res(k)=u1
            enddo

C come back here if we failed to place this conformation - note that we
C had to pull this statement out of the inot_made_yet if statment...

809         r=minval(rand_res(1:num_peps(i1,i2)))

C now we need to determine which dipeptide pdb file to read - we could
C either be reading from a library or we could be reading from pre-built
C do the following if reading from the dipeptide library...

            if(inot_made_yet(j).eq.1) then

C make sure that we set i_am_original for this residue to zero...

              i_am_original(j)=0

c             if(nfails_res.eq.0)
c    &        write(*,*)'reading library for res# ',j-1,j,dipep_cur

C v3.0 randomly select whether to put in a helical residue here
C basic form is correct here but we need to be sure my_PSIPRED_vl2 
C is accessed correctly - jkl should be sequence# corresponding to this
C fragment, 1 should always be correct (we only allow one pdb to act as
C the template), and jlock should be the chain-local res# - so it all
C looks correct at this stage but obviously needs testing...

C v3.4 instead we will automatically select helix if the confidence
C value is above the cutoff specified in the input file:
C   helix_cut

C v3.5 implement same idea for sheets...

C note that we don't allow for dipeptides to have HE or EH conformations
C (i.e. helix followed immediately by a sheet) - so, as written, the if
C statement below will give preference to helices: a HE segment will be
C treated as a HC segment, and a EH segment will become a CH segment

Cv3.0         call random_number(v1)
              jkl=my_frags_seq_num(my_frag)

Cv3.0         if(v1.lt.my_PSIPRED_vl2(jkl,1,jloc)) then
Cv3.0           k=1
Cv3.0         else
Cv3.0           k=minloc(rand_res(1:num_peps(i1,i2)),1)
Cv3.0           rand_res(k)=2.0 ! make sure k is not selected again
Cv3.0         endif

              if(my_PSIPRED_vl2(jkl,1,jloc-1).lt.helix_cut.and.
     &           my_PSIPRED_vl2(jkl,1,jloc  ).lt.helix_cut.and.
     &           my_PSIPRED_vl3(jkl,1,jloc-1).lt.sheet_cut.and.
     &           my_PSIPRED_vl3(jkl,1,jloc  ).lt.sheet_cut) then
                k=minloc(rand_res(1:dipep_both_coils_end(i1,i2)),1)
              elseif(my_PSIPRED_vl2(jkl,1,jloc-1).ge.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).lt.helix_cut) then
                k=minloc(rand_res(dipep_helix_Cter_beg(i1,i2):
     &                            dipep_helix_Cter_end(i1,i2)),1)
                k=k+dipep_helix_Cter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl2(jkl,1,jloc-1).ge.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).ge.helix_cut) then
                k=minloc(rand_res(dipep_both_helix_beg(i1,i2):
     &                            dipep_both_helix_end(i1,i2)),1)
                k=k+dipep_both_helix_beg(i1,i2)-1
              elseif(my_PSIPRED_vl2(jkl,1,jloc-1).lt.helix_cut.and.
     &               my_PSIPRED_vl2(jkl,1,jloc  ).ge.helix_cut) then
                k=minloc(rand_res(dipep_helix_Nter_beg(i1,i2):
     &                            dipep_helix_Nter_end(i1,i2)),1)
                k=k+dipep_helix_Nter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc-1).ge.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).lt.sheet_cut) then
                k=minloc(rand_res(dipep_sheet_Cter_beg(i1,i2):
     &                            dipep_sheet_Cter_end(i1,i2)),1)
                k=k+dipep_sheet_Cter_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc-1).ge.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).ge.sheet_cut) then
                k=minloc(rand_res(dipep_both_sheet_beg(i1,i2):
     &                            dipep_both_sheet_end(i1,i2)),1)
                k=k+dipep_both_sheet_beg(i1,i2)-1
              elseif(my_PSIPRED_vl3(jkl,1,jloc-1).lt.sheet_cut.and.
     &               my_PSIPRED_vl3(jkl,1,jloc  ).ge.sheet_cut) then
                k=minloc(rand_res(dipep_sheet_Nter_beg(i1,i2):
     &                            dipep_sheet_Nter_end(i1,i2)),1)
                k=k+dipep_sheet_Nter_beg(i1,i2)-1
              endif

              rand_res(k)=2.0 ! make sure k is not selected again

              xm=0.0 ! remember to zero out all entries
              ym=0.0 ! remember to zero out all entries
              zm=0.0 ! remember to zero out all entries
              xm(1)=x_dipep(dipep_N1(i1,i2),k,i1,i2)
              ym(1)=y_dipep(dipep_N1(i1,i2),k,i1,i2)
              zm(1)=z_dipep(dipep_N1(i1,i2),k,i1,i2)
              xm(2)=x_dipep(dipep_CA1(i1,i2),k,i1,i2)
              ym(2)=y_dipep(dipep_CA1(i1,i2),k,i1,i2)
              zm(2)=z_dipep(dipep_CA1(i1,i2),k,i1,i2)
              xm(3)=x_dipep(dipep_C1(i1,i2),k,i1,i2)
              ym(3)=y_dipep(dipep_C1(i1,i2),k,i1,i2)
              zm(3)=z_dipep(dipep_C1(i1,i2),k,i1,i2)
              xm(4)=x_dipep(dipep_O1(i1,i2),k,i1,i2)
              ym(4)=y_dipep(dipep_O1(i1,i2),k,i1,i2)
              zm(4)=z_dipep(dipep_O1(i1,i2),k,i1,i2)

c      write(*,*)'BEF1',xm(1),ym(1),zm(1),rf(1,1),rf(1,2),rf(1,3)
c      write(*,*)'BEF2',xm(2),ym(2),zm(2),rf(2,1),rf(2,2),rf(2,3)
c      write(*,*)'BEF3',xm(3),ym(3),zm(3),rf(3,1),rf(3,2),rf(3,3)
c      write(*,*)'BEF4',xm(4),ym(4),zm(4),rf(4,1),rf(4,2),rf(4,3)

              lk=0
              do l=dipep_natom1(i1,i2)+1,dipep_natoms(i1,i2)
                lk=lk+1
                xm(4+lk)=x_dipep(l,k,i1,i2)
                ym(4+lk)=y_dipep(l,k,i1,i2)
                zm(4+lk)=z_dipep(l,k,i1,i2)
              enddo
              nm=4+dipep_natom2(i1,i2)
              itoadd=dipep_natom2(i1,i2) ! store how many atoms we'll add
              itoskp=dipep_natom1(i1,i2) ! store how many atoms to skip

            elseif(inot_made_yet(j).eq.0) then

C do the following if reading the dipeptide from pre-built structure
C note that we have a simple way to switch back to use the dipeptide
C library if we fail excessively with the pre-built structure: we just
C reset inot_made_yet to 1 for this particular dipeptide...

c             if(nfails_res.eq.0)
c    &        write(*,*)'reading strctre for res# ',j-1,j,dipep_cur

              xm=0.0 ! remember to zero out all entries
              ym=0.0 ! remember to zero out all entries
              zm=0.0 ! remember to zero out all entries
              xm(1)=x_temp4(temp4_N(jx-1),jx-1)
              ym(1)=y_temp4(temp4_N(jx-1),jx-1)
              zm(1)=z_temp4(temp4_N(jx-1),jx-1)
              xm(2)=x_temp4(temp4_CA(jx-1),jx-1)
              ym(2)=y_temp4(temp4_CA(jx-1),jx-1)
              zm(2)=z_temp4(temp4_CA(jx-1),jx-1)
              xm(3)=x_temp4(temp4_C(jx-1),jx-1)
              ym(3)=y_temp4(temp4_C(jx-1),jx-1)
              zm(3)=z_temp4(temp4_C(jx-1),jx-1)
              xm(4)=x_temp4(temp4_O(jx-1),jx-1)
              ym(4)=y_temp4(temp4_O(jx-1),jx-1)
              zm(4)=z_temp4(temp4_O(jx-1),jx-1)
c             do l=1,n_temp4(jx-1)
c               xm(4+l)=x_temp4(l,jx-1)
c               ym(4+l)=y_temp4(l,jx-1)
c               zm(4+l)=z_temp4(l,jx-1)
c             enddo
c             nm=4+n_temp4(jx-1)
              do l=1,n_temp4(jx)
                xm(4+l)=x_temp4(l,jx)
                ym(4+l)=y_temp4(l,jx)
                zm(4+l)=z_temp4(l,jx)
              enddo
c             nm=4+n_temp4(jx-1)+n_temp4(jx)
              nm=4+n_temp4(jx)
              itoadd=n_temp4(jx) ! store how many atoms we will add
              itoskp=0           ! store how many atoms to skip

            endif

C at this point we should have our files 'target' and 'temp4' and be
C ready to do a superposition...

            call pdbsup(nm,xm,ym,zm,rf,rms)
c           write(*,*)'rms = ',rms,' nm = ',nm
c           write(91,'(f10.5)')rms

c      write(*,*)'AFT1',xm(1),ym(1),zm(1),rf(1,1),rf(1,2),rf(1,3)
c      write(*,*)'AFT2',xm(2),ym(2),zm(2),rf(2,1),rf(2,2),rf(2,3)
c      write(*,*)'AFT3',xm(3),ym(3),zm(3),rf(3,1),rf(3,2),rf(3,3)
c      write(*,*)'AFT4',xm(4),ym(4),zm(4),rf(4,1),rf(4,2),rf(4,3)

C now decide whether to accept this on the basis of its rms...

            if(rms.gt.rms_crit) then

              nfails_res=nfails_res+1
ccc           write(*,*)'HELP1C!! ',nfails_res

C note that the last criterion ensures that we immediately quit if a
C structured peptide was used to add the new residue

              if(nfails_res.ge.num_peps(i1,i2).or.
     &           nfails_res.ge.nfails_res_max.or.
     &           inot_made_yet(j).eq.0) then

                nfails_chn=nfails_chn+1

C 2023_v1.1 - only extend if nfails_chn_max>0
 
                if(nfails_chn.ge.nfails_chn_max.and.
     &             nfails_chn_max.gt.0) then
                  write(*,*)'rms   extending Cter chain # ',i,rms
                  isentback=1
                  goto 404 ! extend the chain
                else
c                 write(*,*)'rms   restartng Cter chain # ',i,rms

C v4.3 make sure we try to backtrack

                  if(nbacktrack.gt.0) itrybacktrack=1

                  isentback=1 ! BUGFIX ADDED LATE TO SPEED UP
                  goto 505 ! restart the chain
                endif

              else
                goto 809 ! pick another conformation
              endif

            endif

C if we got here then we may be about to include the new amino acid...
C first, get the atoms that we want from the superimpose output

            do m=1,itoadd
              xm(m)=xm(m+4)
              ym(m)=ym(m+4)
              zm(m)=zm(m+4)
              iresm(m)=jx ! v2.1 store res# for debugging
            enddo

C if membrane_dive then get the z coords of the N atoms

            if(membrane_dive) then

C TEMP TEMP TEMP v2.6 may need to switch z1 and z2...
              z1=zm(1+4) ! store z coord of N atom in new res
     &                   ! this is always first atom
              z2=rf(1,3) ! store z coord of prev N atom - 1st fixd atom
ccc           z2=zm(1+4) ! store z coord of N atom in prev res
ccc  &                   ! this is always first atom
ccc           z1=zm(1+4+n_temp4(jx-1)) ! store z coord of 2nd N atom
            endif

C if membrane_dive and this is the first residue we're adding then we
C can use the previous residue to tell us which zone we're in - we've
C got four possibilities:

C (1) i_go_up=1 : tail starts below membrane (must stay that side)
C (2) i_go_up=2 : tail starts within membrane (and must move in neg)
C (3) i_go_up=3 : tail starts within membrane (and must move in pos)
C (4) i_go_up=4 : tail starts above membrane (must stay that side)

            if(membrane_dive.and.j.eq.ibeg_tot_cters(i)) then
              if(zmem_neg2.eq.0.000) then
                if(z2.le.zmem_neg1) then
                  i_go_up=1
                elseif(z2.le.0.0.and.z2.gt.zmem_neg1) then
                  i_go_up=2
                elseif(z2.le.zmem_pos1.and.z2.gt.0.0) then
                  i_go_up=3
                else
                  i_go_up=4
                endif
              else
                if(z2.le.zmem_neg1) then
                  i_go_up=1
                elseif(z2.le.zmem_ave1.and.z2.gt.zmem_neg1) then
                  i_go_up=2
                elseif(z2.le.zmem_pos1.and.z2.gt.zmem_ave1) then
                  i_go_up=3
                elseif(z2.le.zmem_neg2.and.z2.gt.zmem_pos1) then
                  i_go_up=4
                elseif(z2.le.zmem_ave2.and.z2.gt.zmem_neg2) then
                  i_go_up=5
                elseif(z2.le.zmem_pos2.and.z2.gt.zmem_ave2) then
                  i_go_up=6
                else
                  i_go_up=7
                endif
              endif
            endif
            
C now, if membrane_dive, we can ask if we're going in right direction
C note that imem_quit - which is used for quitting this tail conf - is
C only set non-zero  if membrane_dive is true... remember that z1 is the
C new residue and z2 is the old residue...

            imem_quit=0 
            if(membrane_dive) then
              if(zmem_neg2.eq.0.000) then
                if(i_go_up.eq.1) then
                  if(z1.gt.zmem_neg1) imem_quit=1
                elseif(i_go_up.eq.2) then
                  if(z2.ge.zmem_neg1.and.z1.gt.z2) imem_quit=1
                elseif(i_go_up.eq.3) then
                  if(z2.le.zmem_pos1.and.z1.lt.z2) imem_quit=1
                elseif(i_go_up.eq.4) then
                  if(z1.lt.zmem_pos1) imem_quit=1
                endif
              else
                if(i_go_up.eq.1) then
                  if(z1.gt.zmem_neg1) imem_quit=1
                elseif(i_go_up.eq.2) then
                if(z2.ge.zmem_neg1.and.z1.gt.z2) imem_quit=1
              elseif(i_go_up.eq.3) then
                if(z2.le.zmem_pos1.and.z1.lt.z2) imem_quit=1
              elseif(i_go_up.eq.4) then
                if(z1.lt.zmem_pos1) imem_quit=1
                if(z1.gt.zmem_neg2) imem_quit=1
              elseif(i_go_up.eq.5) then
                if(z2.ge.zmem_neg2.and.z1.gt.z2) imem_quit=1
              elseif(i_go_up.eq.6) then
                if(z2.le.zmem_pos2.and.z1.lt.z2) imem_quit=1
              elseif(i_go_up.eq.7) then
                if(z1.lt.zmem_pos2) imem_quit=1
              endif
            endif
c           write(*,*)
c           write(*,*)'check z1 z2 Cter settings ',z1,z2,i_go_up
c           write(*,*)
            endif
   
C now, before finally adding this residue we need to do a clash check
C we'll use the same code to implement the quit due to membrane_dive
C rather than adding in all the same code above

C TEMP TEMP TEMP we need to remember that we are currently only checking
C against fixed residues - this will be everything except HETATM entries
C and except for other residues in the tail....

            iclash=0
            do m1=1,itoadd
              xt=xm(m1)
              yt=ym(m1)
              zt=zm(m1)

C skip clash-checking if this atom is not on the grid - note that this
C used to be a fatal error when I made a smaller clash grid that only
C encompassed the atoms near to the beginning of the tail...

              if(xt.lt.xmin.or.xt.gt.xmax.or.
     &           yt.lt.ymin.or.yt.gt.ymax.or.
     &           zt.lt.zmin.or.zt.gt.zmax) cycle

              it1=int((xt-xmin)*xinv)+1
              jt1=int((yt-ymin)*yinv)+1
              kt1=int((zt-zmin)*zinv)+1
              do m2=1,temp_grid(it1,jt1,kt1)

C skip the grid atom if it's the jx-1'th residue of same chain

                if(cchn_grid(m2,it1,jt1,kt1).eq.
     &             cur_chainID(my_frag).and.
     &             rchn_grid(m2,it1,jt1,kt1).eq.jx-1) then
c                 write(*,*)'skip clash with ',rchn_grid(m2,it1,jt1,kt1),
c    &                                         cchn_grid(m2,it1,jt1,kt1)
                  cycle
                endif

C otherwise, measure the distance and look for a clash

                xu=xatm_grid(m2,it1,jt1,kt1)
                yu=yatm_grid(m2,it1,jt1,kt1)
                zu=zatm_grid(m2,it1,jt1,kt1)
                dist2=((xt-xu)**2+(yt-yu)**2+(zt-zu)**2)

C v4.2 do this if it's a potential ATOM-ATOM clash 

                if(ihet_grid(m2,it1,jt1,kt1).eq.0) then
                  if(dist2.lt.clash_aa2*clash_fac2) then
                    iclash=1
                    jclash=rchn_grid(m2,it1,jt1,kt1)
c                   write(*,*)'intermolec clash ',iresm(m1),
c    &                rchn_grid(m2,it1,jt1,kt1),sqrt(dist2)
                    goto 9151
                  endif

C v4.2 do this if it's a potential ATOM-HETATM clash 

                elseif(ihet_grid(m2,it1,jt1,kt1).eq.1) then
                  if(dist2.lt.clash_ha2*clash_fac2) then
                    iclash=1
                    jclash=rchn_grid(m2,it1,jt1,kt1)
c                   write(*,*)'intermolec clash ',iresm(m1),
c    &                rchn_grid(m2,it1,jt1,kt1),sqrt(dist2)
                    goto 9151
                  endif
                endif

              enddo
            enddo

C check for intra-tail clashes - note that itochk contains residues j-1
C downwards so we need to explicitly check to ignore j-1:j clashes
C note also that it's probably quicker to look for clashes starting with
C the last atoms in the "itochk" list - these will be the most likely to
C be clashing with the current residue, and if we're going to fail we
C want to fail quickly...

            do m1=1,itoadd
              xt=xm(m1)
              yt=ym(m1)
              zt=zm(m1)
              do m2=itochk,1,-1 ! v2.7 see comment above
                xu=xtail(m2)
                yu=ytail(m2)
                zu=ztail(m2)
                dist2=((xt-xu)**2+(yt-yu)**2+(zt-zu)**2)
C v4.2          if(dist2.lt.clash2) then
                if(dist2.lt.clash_aa2*clash_fac2) then
                  iclash=1

C v2.7 here we ignore clashes between adjacent residues

                  jclash=ires_tot_tail(m2)
                  if(abs(jclash-j).eq.1) then
                    iclash=0
                    cycle
                  endif

C we have to be careful not to quit if the apparent "clash" is between
C two residues that are part of the same structured domain...
C here, remember that "j" and "jclash" are both in total-residue count
C numbering so we can plug them directly in and find their domains -
C note also that they must, by definition, be part of the same molecule

C v4.4 - we also need to check that they still have the inot_made_yet
C status set to 0 - these could have been updated to 1 if we failed to
C build the tail many times... - add this extra condition to the if
C statement so we only come here if all these things are true:

c                 if(i_in_dom_fasta_tot(j).eq.
c    &               i_in_dom_fasta_tot(jclash).and.
c    &               i_in_dom_fasta_tot(j).gt.0) then

                  if(i_in_dom_fasta_tot(j).eq.
     &               i_in_dom_fasta_tot(jclash).and.
     &               i_in_dom_fasta_tot(j).gt.0.and.
     &               inot_made_yet(j).eq.0.and.
     &               inot_made_yet(jclash).eq.0) then

                    iclash=0 ! set back to zero to ignore the clash
                    cycle
c                   write(*,*)'ignoring intra-tail clash ',
c    &                         iresm(m1),iresn(m2),sqrt(dist2),
c    &                         i_in_dom_fasta_tot(j)

C v4.5 we must do something similar for ATOM-HETATM entries where the
C HETATM is associated with the tail... TEMP TEMP TEMP

                  elseif(i_in_dom_fasta_tot(j).eq.dtail(m2).and.
     &                   inot_made_yet(j).eq.0) then
                    iclash=0
                    cycle

                  else
c                   write(*,*)'deadly intra-tail clash ',
c    &                         iresm(m1),iresn(m2),sqrt(dist2)
c                   write(*,*)'xm,ym,zm = ',xt,yt,zt
c                   write(*,*)'xn,yn,zn = ',xu,yu,zu
                    goto 9151
                  endif
                endif
              enddo
            enddo

9151        continue ! come here if we had a clash...

            if(iclash.eq.1.or.imem_quit.eq.1) then

C v4.7 6/10/2022 if imem_quit=1 let's restart the chain no questions
C asked - if we don't do this then we might keep snaking along the
C membrane leading to a highly biased distribution...
C UPDATE - the code seems to hang if we do this :(

ccc           if(imem_quit.eq.1) then
ccc             goto 5432   ! restart the chain
ccc           endif

              nfails_res=nfails_res+1
ccc           write(*,*)'HELP2C!! ',nfails_res

              if(nfails_res.ge.num_peps(i1,i2).or.
     &           nfails_res.ge.nfails_res_max.or.
     &           inot_made_yet(j).eq.0) then

                nfails_chn=nfails_chn+1

C 2023_v1.1 - only extend if nfails_chn_max>0
 
                if(nfails_chn.ge.nfails_chn_max.and.
     &             nfails_chn_max.gt.0) then
                  write(*,*)'clash extending Cter chain # ',i,j,jclash,
     &                       sqrt(dist2)
                  isentback=1
                  goto 404 ! extend the chain
                else
c                 write(*,*)'clash restartng Cter chain # ',i,j,jclash,
c    &                       sqrt(dist2),nfails_chn

C v4.3 make sure we try to backtrack

                  if(nbacktrack.gt.0) itrybacktrack=1

                  isentback=1 ! BUGFIX ADDED LATE TO SPEED UP
                  goto 505 ! restart the chain
                endif

              else
                goto 809 ! pick another conformation
              endif

            endif ! done with the "if clashing..."

C if we got here then we really are adding this residue...
C now add the atoms of the new residue to string66_tail...

C first, add on the other residues from current_tail.pdb...

            do ito=1,itochk
              string66_tail_tmp(ito)=string66_tail(ito)
              ires_tot_tail_tmp(ito)=ires_tot_tail(ito)
              dtail_tmp(ito)=dtail(ito)
              xtail_tmp(ito)=xtail(ito)
              ytail_tmp(ito)=ytail(ito)
              ztail_tmp(ito)=ztail(ito)
            enddo

C now, write out the new residue (it's more C-terminal than others)
C note that we write it out with the correct res# and chainID
C note also that we store the N,CA,C,O coordinates...

C note also that if we're using a dipeptide we need to make sure that
C we're looking only in the *second* of the two amino acids for the atom
C names - hence the reason we add "itoskp" to l for string66_dipep...

            if(inot_made_yet(j).eq.1) then
              do l=1,itoadd
                m=l+itoskp
                if(string66_dipep(m,k,i1,i2)(13:16).eq.' N  ') then
                  rf_prev(1,1)=xm(l)
                  rf_prev(1,2)=ym(l)
                  rf_prev(1,3)=zm(l)
                elseif(string66_dipep(m,k,i1,i2)(13:16).eq.' CA ') then
                  rf_prev(2,1)=xm(l)
                  rf_prev(2,2)=ym(l)
                  rf_prev(2,3)=zm(l)
                elseif(string66_dipep(m,k,i1,i2)(13:16).eq.' C  ') then
                  rf_prev(3,1)=xm(l)
                  rf_prev(3,2)=ym(l)
                  rf_prev(3,3)=zm(l)
                elseif(string66_dipep(m,k,i1,i2)(13:16).eq.' O  ') then
                  rf_prev(4,1)=xm(l)
                  rf_prev(4,2)=ym(l)
                  rf_prev(4,3)=zm(l)
                endif

C v2.7 store information in internal strings to avoid constant writing

                write(string66_tail_tmp(l+itochk),81)
     &                      string66_dipep(m,k,i1,i2)(1:21),
     &                      cur_chainID(my_frag),
     &                      jx,xm(l),ym(l),zm(l),
     &                      string66_dipep(m,k,i1,i2)(55:66)     ! v2.7
                ires_tot_tail_tmp(l+itochk)=j
                dtail_tmp(l+itochk)=dm(l)
                xtail_tmp(l+itochk)=xm(l)
                ytail_tmp(l+itochk)=ym(l)
                ztail_tmp(l+itochk)=zm(l)

              enddo
            else
              do l=1,itoadd

C v4.5 since each residue might now have associated with it a bunch of
C HETATM entries we need to make sure that we're only looking at ATOMs

                if(string66_temp4(l,jx)(1:4).eq.'ATOM') then

                  if(string66_temp4(l,jx)(13:16).eq.' N  ') then
                    rf_prev(1,1)=xm(l)
                    rf_prev(1,2)=ym(l)
                    rf_prev(1,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' CA ') then
                    rf_prev(2,1)=xm(l)
                    rf_prev(2,2)=ym(l)
                    rf_prev(2,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' C  ') then
                    rf_prev(3,1)=xm(l)
                    rf_prev(3,2)=ym(l)
                    rf_prev(3,3)=zm(l)
                  elseif(string66_temp4(l,jx)(13:16).eq.' O  ') then
                    rf_prev(4,1)=xm(l)
                    rf_prev(4,2)=ym(l)
                    rf_prev(4,3)=zm(l)
                  endif

                endif

C v2.7 store information in internal strings to avoid constant writing
C v4.5 make sure we don't change the residue number of HETATMs

                if(string66_temp4(l,jx)(1:4).eq.'ATOM') then

                  write(string66_tail_tmp(l+itochk),81)
     &                      string66_temp4(l,jx)(1:21),       
     &                      cur_chainID(my_frag),
     &                      jx,xm(l),ym(l),zm(l),
     &                      string66_temp4(l,jx)(55:66)           ! v2.7

                elseif(string66_temp4(l,jx)(1:4).ne.'ATOM') then

                  write(string66_tail_tmp(l+itochk),94)
     &                      string66_temp4(l,jx)(1:21),       
ccc  &                      cur_chainID(my_frag), ! v4.6 make sure
     &                      fin_chainID(my_frag), ! chainID is final
     &                      string66_temp4(l,jx)(23:26),     
     &                      xm(l),ym(l),zm(l),
     &                      string66_temp4(l,jx)(55:66)           ! v2.7

                endif

                ires_tot_tail_tmp(l+itochk)=j
                dtail_tmp(l+itochk)=dm(l)
                xtail_tmp(l+itochk)=xm(l)
                ytail_tmp(l+itochk)=ym(l)
                ztail_tmp(l+itochk)=zm(l)

              enddo

            endif

C v2.7 here we should put the new atoms on the grid for clash checking




C update itochk to include the new residue

            itochk=itochk+itoadd

C finally, copy back to string66_tail...

            do ito=1,itochk
              string66_tail(ito)=string66_tail_tmp(ito)
              ires_tot_tail(ito)=ires_tot_tail_tmp(ito)
              dtail(ito)=dtail_tmp(ito)
              xtail(ito)=xtail_tmp(ito)
              ytail(ito)=ytail_tmp(ito)
              ztail(ito)=ztail_tmp(ito)
            enddo

C if this is the final residue then we write out to current_tail.pdb

            if(j.eq.iend_tot_cters(i)) then

              ihaveoxt=0
              open(unit=12,file='current_tail.pdb',status='unknown')
              do l=1,itochk
                write(12,'(a66)')string66_tail(l)
                if(string66_tail(l)(14:16).eq.'OXT') ihaveoxt=1
              enddo

C and we add the OXT atom here as the last atom in the file
C v4.5 we should only do this if we don't already have an OXT atom in
C the file - see the check above - OXT atoms seem to be possible
C additions at multiple points in the code :(

              if(ihaveoxt.eq.0) then

              write(*,*)'WRITING OXT ATOM TO CTER TAIL'

              a(1)=rf_prev(4,1)
              a(2)=rf_prev(4,2)
              a(3)=rf_prev(4,3)
              b(1)=rf_prev(2,1)
              b(2)=rf_prev(2,2)
              b(3)=rf_prev(2,3)
              c(1)=rf_prev(3,1)
              c(2)=rf_prev(3,2)
              c(3)=rf_prev(3,3)

              bondcur=1.25
              anglcur=118.0*3.14159265/180.0
              dihecur=180.0*3.14159265/180.0
              st=sin(anglcur)
              zpd=-bondcur*cos(anglcur)
              xpd= bondcur*cos(dihecur)*st
              ypd= bondcur*sin(dihecur)*st
              do l=1,3
                vca(l)=a(l)-c(l)
                vcb(l)=b(l)-c(l)
              enddo
              call acrosb(vca,vcb,yp)
              call acrosb(vcb,yp,xp)
              call acrosb(xp,yp,zp)
              xoxt=xp(1)*xpd+yp(1)*ypd+zp(1)*zpd+c(1)
              yoxt=xp(2)*xpd+yp(2)*ypd+zp(2)*zpd+c(2)
              zoxt=xp(3)*xpd+yp(3)*ypd+zp(3)*zpd+c(3)
              write(12,89)string66_tail(itochk)(1:12),
     &                    string66_tail(itochk)(17:21),
     &                    cur_chainID(my_frag),
     &                    string66_tail(itochk)(23:26),
     &                    xoxt,yoxt,zoxt,string66_tail(itochk)(55:66)
89            format(a12,' OXT',a5,a,a4,4x,3f8.3,a12)

              elseif(ihaveoxt.eq.1) then

              write(*,*)'NO NEED TO ADD OXT ATOM TO CTER TAIL'

              endif

C finally, close current_tail.pdb

              close(12)

            endif

          enddo

        endif

C at this point we should have loopy_inp.pdb and a completed 
C current_tail.pdb - we want to add the tail coords found in 
C current_tail.pdb into loopy_inp.pdb and then copy back to 
C loopy_inp.pdb - so we'll first make a temporary copy:

C v4.5 I don't trust the system copy commands at all so omit them

ccc     call system('cp loopy_inp.pdb loopy_tmp.pdb')
ccc     write(*,*)'istatus 001 = ',istatus
        open(unit=11,file='loopy_inp.pdb',status='unknown') ! input
        open(unit=12,file='loopy_tmp.pdb',status='unknown') ! output
8501    read(11,'(a80)',end=8502)char80
        write(12,'(a80)')char80
        goto 8501
8502    close(11)
        close(12)

C v3.1 keep a copy of the tail

        write(char4,'(i4)')my_frag
        if(char4(1:1).eq.' ') char4(1:1)='0'
        if(char4(2:2).eq.' ') char4(2:2)='0'
        if(char4(3:3).eq.' ') char4(3:3)='0'
        if(char4(4:4).eq.' ') char4(4:4)='0'

        if(typ_tot_tails(itmp).eq.1) then

          tailname(1:26)='current_tail_Nter_XXXX.pdb'
          tailname(19:22)=char4
          open(unit=12,file='current_tail.pdb',status='unknown')
          open(unit=72,file=tailname,status='unknown')
4881      read(12,'(a80)',end=4882)char80
          write(72,'(a80)')char80
          goto 4881
4882      close(12)
          close(72)
          write(*,*)
 
        elseif(typ_tot_tails(itmp).eq.2) then

          tailname(1:26)='current_tail_Cter_XXXX.pdb'
          tailname(19:22)=char4
          open(unit=12,file='current_tail.pdb',status='unknown')
          open(unit=72,file=tailname,status='unknown')
4883      read(12,'(a80)',end=4884)char80
          write(72,'(a80)')char80
          goto 4883
4884      close(12)
          close(72)
          write(*,*)

        endif
 
        open(unit=11,file='loopy_tmp.pdb',status='unknown') ! input
        open(unit=12,file='current_tail.pdb',status='unknown') ! input
        open(unit=13,file='loopy_inp.pdb',status='unknown') ! output

C this will assume that we never have to replace this tail which
C might hurt us in some cases...

        if(typ_tot_tails(itmp).eq.1) then

C get the correct values of jbeg and jend at this point

          jbeg=ibeg_loc_nters(i)+cur_offset(my_frag)
          jend=iend_loc_nters(i)+cur_offset(my_frag)

          iatpoint=0
9027      read(11,'(a66)',end=9028)string66
          read(string66,'(21x,a,i4)')chainID,jkl

C make sure that we skip residues present in loopy_tmp.pdb that are
C going to be replaced by those in the tail just built
C v4.5 - do this only if it's an ATOM

          if(chainID.eq.ichn_tot_nters(i).and.
     &       string66(1:4).eq.'ATOM'.and. ! v4.5
     &       jkl.ge.jbeg.and.jkl.le.jend) goto 9027

C otherwise, read until we get to the first atom of the res *after* the tail
C v4.5 - do this only if it's an ATOM

          if(chainID.eq.ichn_tot_nters(i).and.
     &       string66(1:4).eq.'ATOM'.and.jkl.eq.jend+1) then
            if(iatpoint.eq.0) then
9029          read(12,'(a66)',end=9030)strong66
              if(strong66(1:4).ne.'ATOM') goto 9029 ! v4.5 skip HETATMs
              read(strong66,'(21x,a,i4)')chainID,jkl
              if(chainID.eq.ichn_tot_nters(i).and.jkl.ge.jbeg.and.
     &           jkl.le.jend) write(13,'(a66)')strong66
              goto 9029
9030          rewind(12) ! v4.5 we used to close, now we just rewind
              write(13,'(a66)')string66
              iatpoint=1
            else
              write(13,'(a66)')string66
            endif
          else
            write(13,'(a66)')string66
          endif
          goto 9027
9028      close(11)

C v4.5 at this point we can add on any tail-associated HETATM entries in
C the current_tail.pdb...

9025      read(12,'(a66)',end=9026)strong66
          if(strong66(1:4).ne.'HETA') goto 9025 ! v4.5 skip ATOMs
          write(13,'(a66)')strong66
          goto 9025
9026      close(12)

          close(13)
ccc       call system('rm current_tail.pdb')

          write(*,*)'added Nter chain # ',i,' time to get excited?'

        elseif(typ_tot_tails(itmp).eq.2) then

C get the correct values of jbeg and jend at this point

          jbeg=ibeg_loc_cters(i)+cur_offset(my_frag)
          jend=iend_loc_cters(i)+cur_offset(my_frag)

C get the chainID type of the current chain - note that this is needed
C for Cter tails but not Nter tails... (comment added in v3.8+)

          do nj=1,chainIDmax
            if(ichn_tot_cters(i).eq.blah(nj:nj)) then
              nnntmp=nj
              exit
            endif
          enddo

          iatpoint=0
9127      read(11,'(a66)',end=9128)string66
          read(string66,'(21x,a,i4)')chainID,jkl

C make sure that we skip residues present in loopy_tmp.pdb that are
C going to be replaced by those in the tail just built
C v4.5 - do this only if it's an ATOM

          if(chainID.eq.ichn_tot_cters(i).and.
     &       string66(1:4).eq.'ATOM'.and. ! v4.5
     &       jkl.ge.jbeg.and.jkl.le.jend) goto 9127

C otherwise, read until we get to the first atom of the res *after* the
C tail - note that the first line of the if statement is
C straightforward, the second line is there to catch cases where the 
C next residue has a different chainID - we handle this by comparing the 
C chain type of current tail (nnntmp) with the chain type of the atom
C that was just read from loopy_tmp.pdb (mmmtmp) - if the latter is one
C greater than the former then we can also invoke the writing of the
C current tail residue...
C v4.5 - do this only if it's an ATOM
C v4.5 - change the requirement from "mmmtmp.eq.nnntmp+1"
C        to "mmmtmp.ne.nnntmp" - doing this for L12 tails where
C        the added tail is chain A and next (fake GLY) res is chain E
C        - the pre-existing code doesn't catch it cos E<>B...

C get the chainIDtype of the current line from loopy_tmp.pdb

          do nj=1,chainIDmax
            if(chainID.eq.blah(nj:nj)) then
              mmmtmp=nj
              exit
            endif
          enddo

          if((chainID.eq.ichn_tot_cters(i).and.
     &       string66(1:4).eq.'ATOM'.and.jkl.gt.jend).or.
     &       (mmmtmp.ne.nnntmp.and.string66(1:4).eq.'ATOM')) then
ccc  &       (mmmtmp.eq.nnntmp+1.and.string66(1:4).eq.'ATOM')) then
            if(iatpoint.eq.0) then
9129          read(12,'(a66)',end=9130)strong66
              if(strong66(1:4).ne.'ATOM') goto 9129 ! v4.5 skip HETATMs
              read(strong66,'(21x,a,i4)')chainID,jkl
              if(chainID.eq.ichn_tot_cters(i).and.jkl.ge.jbeg.and.
     &           jkl.le.jend) write(13,'(a66)')strong66
              goto 9129
9130          rewind(12) ! v4.5 we used to close, now we just rewind
              write(13,'(a66)')string66
              iatpoint=1
            else
              write(13,'(a66)')string66
            endif
          else
            write(13,'(a66)')string66
          endif
          goto 9127

9128      close(11)

C v4.5 we need to make sure that the C-terminal tail is added on
C correctly when it's the very last part of the protein - we won't have
C found it with the code above since that's looking for the first atom
C *after* the tail in loopy_tmp.pdb so we'll do a final check here and
C make sure to add on the tail if we haven't already done it

          if(iatpoint.eq.0) then
            write(*,*)'this is the last C-term tail added here'
9229        read(12,'(a66)',end=9230)strong66
            read(strong66,'(21x,a,i4)')chainID,jkl
            if(chainID.eq.ichn_tot_cters(i).and.jkl.ge.jbeg.and.
     &         jkl.le.jend) write(13,'(a66)')strong66
            goto 9229
9230        rewind(12) ! don't close as we may still want HETATMs...
            iatpoint=1
          endif

C v4.5 at this point we can add on any tail-associated HETATM entries in
C the current_tail.pdb...

9125      read(12,'(a66)',end=9126)strong66
          if(strong66(1:4).ne.'HETA') goto 9125 ! v4.5 skip ATOMs
          write(13,'(a66)')strong66
          goto 9125
9126      close(12)

          close(13)
ccc       call system('rm current_tail.pdb')

          write(*,*)'added Cter chain # ',i,' time to get excited?'

        endif

C now work on the next tail

      enddo 

C now read the final pdb and put in the temperature factors
C and change the residue numbers so they are intra-chain numbers

C note that sys_string9b will be the CARTOON.pdb...

      write(*,*)
      write(*,*)'almost done: writing final pdb file'
      write(*,*)
      mmm=len(trim(sys_string9))
      mmm2=len(trim(sys_string9b))
      mmm3=len(trim(sys_string9c)) ! v4.6
      nats=0
      nats2=0

      nmodels=nmodels+1

      open(unit=11,file='loopy_inp.pdb',status='unknown')

      if(nmodels.eq.1) then
        open(unit=32,file=sys_string9(1:mmm),status='unknown')
        open(unit=22,file=sys_string9b(1:mmm2),status='unknown')

C v4.6 if desired we can write all fixed (non-moving) HETATMs to a
C separate file *once* instead of writing to AHEMODEL.pdb...

        if(separate_hetatms) then
          open(unit=23,file=sys_string9c(1:mmm3),status='unknown')
        endif

      else
        goto 4321
      endif

      open(unit=25,file='final_sequence_alignments.txt',
     &     status='unknown')

C first(!), write out the user name and directory

      write(32,'("REMARK:")')
      write(32,'("REMARK:   made using  :  ahemodel_2023_v1.1.exe")')
      mo1=len(trim(my_name))
      mo2=len(trim(my_directory))
      write(my_fmt2,'(a,i0,a)'),'("REMARK:   made by user:  ",a',mo1,')'
      write(32,my_fmt2)my_name(1:mo1)
      write(my_fmt2,'(a,i0,a)'),'("REMARK:   in directory:  ",a',mo2,')'
      write(32,my_fmt2)my_directory(1:mo2)
      write(32,'("REMARK:")')
      write(32,'("REMARK:   input_file used follows:")')
      write(32,'("REMARK:")')

C first(!), write out the input file - after do loop we write out the
C third-to-last line twice so we get a nice horizontal line end...

      do ni=1,ninput_lines
        mo1=len(trim(input_line(ni)))
        mo1=max(mo1,1) ! needed with gfortran to avoid a0 format
        write(my_fmt2,'(a,i0,a)'),'("REMARK:   ",a',mo1,')'
        write(32,my_fmt2)input_line(ni)(1:mo1)
      enddo
      ni=ninput_lines-2
      mo1=len(trim(input_line(ni)))
      mo1=max(mo1,1) ! needed with gfortran to avoid a0 format
      write(my_fmt2,'(a,i0,a)'),'("REMARK:   ",a',mo1,')'
      write(32,my_fmt2)input_line(ni)(1:mo1)
      write(32,'("REMARK:")')

C first, write out any prior REMARK_AHE entries...

      do iii=1,nremarks
        write(32,'(a80)')strang(iii)
      enddo

C first, write out the names of the alignments that were used

      nchains_tmp=0
      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            if(ialign(iii,nnn,mmm).ne.1) cycle
            kkk=len(trim(falign(iii,nnn,mmm)))

C in v1.3 we replace the old write statement - which wasn't very
C friendly - with one that explicitly says which chain ID corresponds to
C which uniprot ID - remember that the final AHEMODEL.pdb may have the
C chain IDs reordered relative to the input TEMPLATE.pdb...

ccc         write(my_fmt2,'(a,i0,a)'),'("REMARK_AHE3:",3i8,2x,a',kkk,')'
ccc         write(32,my_fmt2)iii,nnn,mmm,falign(iii,nnn,mmm)(1:kkk)
            nchains_tmp=nchains_tmp+1
            ijk=mod(nchains_tmp,chainIDmax)
            if(ijk.eq.0) ijk=chainIDmax
            write(32,817)blah(ijk:ijk),falign(iii,nnn,mmm)(1:6),
     &                   nchains_tmp
          enddo
        enddo
      enddo
817   format('REMARK_AHE3: final chain ID ',a,' has uniprot ID ',a6,
     &       ' is chain # ',i6)
      write(32,'(a16)')"REMARK_AHE4: END" ! write this to terminate the list

C second, write out the actual alignments that were used...

815   read(25,'(a100)',end=816)string(1:100)
      write(32,'("REMARK ",a100)')string(1:100)
      goto 815
816   close(25)

C v2.7 write out the alignment statistics from fort.99

      open(unit=99,file='alignment_stats.txt',status='unknown')
9448  read(99,'(a114)',end=9449)string(1:114)
      write(32,'(a114)')string(1:114)
      goto 9448
9449  close(99)

C third, write out SEQRES entries for all chains so that we can use the
C final_pdb as a template for another round if necessary...

      nchains_tmp=0
      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            if(ialign(iii,nnn,mmm).ne.1) cycle
            nchains_tmp=nchains_tmp+1

            do l=1,70
              outstring(l:l)=' '
            enddo
            outstring(1:19)='SEQRES   1 X  123  '
            l=mod(nchains_tmp,chainIDmax)
            if(l.eq.0) l=chainIDmax
            outstring(12:12)=blah(l:l)      

C make sure SEQRES line has correct #res 

            write(outstring(13:17),'(i5)')num_res_in_seq(iii)

            nline=0
            do mp=1,num_res_in_seq(iii),13
              nline=nline+1
              write(outstring(7:10),'(i4)')nline
              ibeg=mp
              iend=min(num_res_in_seq(iii),mp+12)
              j=0
              do i=ibeg,iend
                j=j+1
                if(seq_of_fasta(iii)(i:i).eq.'X') rnam='XXX'
                if(seq_of_fasta(iii)(i:i).eq.'A') rnam='ALA'
                if(seq_of_fasta(iii)(i:i).eq.'C') rnam='CYS'
                if(seq_of_fasta(iii)(i:i).eq.'D') rnam='ASP'
                if(seq_of_fasta(iii)(i:i).eq.'E') rnam='GLU'
                if(seq_of_fasta(iii)(i:i).eq.'F') rnam='PHE'
                if(seq_of_fasta(iii)(i:i).eq.'G') rnam='GLY'
                if(seq_of_fasta(iii)(i:i).eq.'H') rnam='HIS'
                if(seq_of_fasta(iii)(i:i).eq.'I') rnam='ILE'
                if(seq_of_fasta(iii)(i:i).eq.'K') rnam='LYS'
                if(seq_of_fasta(iii)(i:i).eq.'L') rnam='LEU'
                if(seq_of_fasta(iii)(i:i).eq.'M') rnam='MET'
                if(seq_of_fasta(iii)(i:i).eq.'N') rnam='ASN'
                if(seq_of_fasta(iii)(i:i).eq.'P') rnam='PRO'
                if(seq_of_fasta(iii)(i:i).eq.'Q') rnam='GLN'
                if(seq_of_fasta(iii)(i:i).eq.'R') rnam='ARG'
                if(seq_of_fasta(iii)(i:i).eq.'S') rnam='SER'
                if(seq_of_fasta(iii)(i:i).eq.'T') rnam='THR'
                if(seq_of_fasta(iii)(i:i).eq.'V') rnam='VAL'
                if(seq_of_fasta(iii)(i:i).eq.'W') rnam='TRP'
                if(seq_of_fasta(iii)(i:i).eq.'Y') rnam='TYR'
                i1=20+(j-1)*4
                i2=22+(j-1)*4
                outstring(i1:i2)=rnam
              enddo

C remember to blank out the rest of the last line

              do j=i2+1,70
                outstring(j:j)=' '
              enddo

              write(32,'(a70)')outstring(1:70)

            enddo
          enddo
        enddo 
      enddo

C now figure out the domain numbers for each fragments - note that the
C code as written *may* not properly treat cases where there were
C *already* rigid domains defined (i.e. cases where we are trying to
C build up an AHEMODEL structure iteratively) - need to check later
C note that my_domain was defined earlier at the same time that
C i_am_original was defined...

      my_domain_max=-1
      do n=1,num_frags
        ntmp=num_domains(n)

C for all residues initially assign my_domain from i_in_dom_fasta_tot
C for all residues made by loopy, set my_domain=0...
C for all residues with i_am_original=0, set my_domain=0...

        do m1=my_frst_res(n),my_last_res(n)
          my_domain(m1)=i_in_dom_fasta_tot(m1)
          if(i_am_by_loopy(m1).eq.1) my_domain(m1)=0
          if(i_am_original(m1).eq.0) my_domain(m1)=0
          my_domain_max=max(my_domain_max,my_domain(m1))

C v2.0 here we only change my_bfactr if, in the template_pdb, it's equal
C to -1.0 - if it's anything else then we assume that it's from a
C previous iteration of ahemodel and so should be retained unchanged...

          if(my_bfactr(m1).eq.-1.0) then
            if(i_am_original(m1).eq.3) my_bfactr(m1)=99.99
            if(i_am_original(m1).eq.2) my_bfactr(m1)=75.00
            if(i_am_original(m1).eq.1) my_bfactr(m1)=50.00
            if(i_am_original(m1).eq.0) my_bfactr(m1)= 0.00
          endif

        enddo

      enddo

C v3.0 write out HELIX and SHEET definitions...

      nchains_tmp=0
      do iii=1,num_seqs
        do nnn=1,num_pdbs
          do mmm=1,num_chns(nnn)
            if(ialign(iii,nnn,mmm).ne.1) cycle
            kkk=len(trim(falign(iii,nnn,mmm)))
            nchains_tmp=nchains_tmp+1
            ijk=mod(nchains_tmp,chainIDmax)
            if(ijk.eq.0) ijk=chainIDmax

            ires_max=0
            do i=1,9999
              dssp_tmp(i:i)=' '
            enddo            
            open(unit=93,file='DSSP_vs_PSIPRED.txt',status='unknown')
7997        read(93,'(4x,i5,6x,i6,10x,a)',end=7998)ikk,ires,char1
            if(ikk.eq.iii) dssp_tmp(ires:ires)=char1
            if(ikk.eq.iii) ires_max=max(ires,ires_max)
            goto 7997
7998        close(93)

C now find the helices first

            nhelix=0
            i_in_helix=0
            do i=1,ires_max
              if(dssp_tmp(i:i).eq.'H') then
                if(i_in_helix.eq.0) then
                  nhelix=nhelix+1
                  i_in_helix=1
                  ibeg=i
                  iend=i
                elseif(i_in_helix.eq.1) then
                  iend=i
                endif
              else
                if(i_in_helix.eq.1) then
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'A') rnam1='ALA'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'C') rnam1='CYS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'D') rnam1='ASP'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'E') rnam1='GLU'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'F') rnam1='PHE'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'G') rnam1='GLY'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'H') rnam1='HIS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'I') rnam1='ILE'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'K') rnam1='LYS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'L') rnam1='LEU'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'M') rnam1='MET'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'N') rnam1='ASN'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'P') rnam1='PRO'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'Q') rnam1='GLN'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'R') rnam1='ARG'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'S') rnam1='SER'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'T') rnam1='THR'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'V') rnam1='VAL'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'W') rnam1='TRP'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'Y') rnam1='TYR'
                if(seq_of_fasta(iii)(iend:iend).eq.'A') rnam2='ALA'
                if(seq_of_fasta(iii)(iend:iend).eq.'C') rnam2='CYS'
                if(seq_of_fasta(iii)(iend:iend).eq.'D') rnam2='ASP'
                if(seq_of_fasta(iii)(iend:iend).eq.'E') rnam2='GLU'
                if(seq_of_fasta(iii)(iend:iend).eq.'F') rnam2='PHE'
                if(seq_of_fasta(iii)(iend:iend).eq.'G') rnam2='GLY'
                if(seq_of_fasta(iii)(iend:iend).eq.'H') rnam2='HIS'
                if(seq_of_fasta(iii)(iend:iend).eq.'I') rnam2='ILE'
                if(seq_of_fasta(iii)(iend:iend).eq.'K') rnam2='LYS'
                if(seq_of_fasta(iii)(iend:iend).eq.'L') rnam2='LEU'
                if(seq_of_fasta(iii)(iend:iend).eq.'M') rnam2='MET'
                if(seq_of_fasta(iii)(iend:iend).eq.'N') rnam2='ASN'
                if(seq_of_fasta(iii)(iend:iend).eq.'P') rnam2='PRO'
                if(seq_of_fasta(iii)(iend:iend).eq.'Q') rnam2='GLN'
                if(seq_of_fasta(iii)(iend:iend).eq.'R') rnam2='ARG'
                if(seq_of_fasta(iii)(iend:iend).eq.'S') rnam2='SER'
                if(seq_of_fasta(iii)(iend:iend).eq.'T') rnam2='THR'
                if(seq_of_fasta(iii)(iend:iend).eq.'V') rnam2='VAL'
                if(seq_of_fasta(iii)(iend:iend).eq.'W') rnam2='TRP'
                if(seq_of_fasta(iii)(iend:iend).eq.'Y') rnam2='TYR'
                write(32,'("HELIX",i5,i4,1x,a3,1x,a,i5,2x,a3,1x,a,i5)')
     &                          nhelix,nhelix,rnam1,blah(ijk:ijk),ibeg,
     &                                        rnam2,blah(ijk:ijk),iend
                write(22,'("HELIX",i5,i4,1x,a3,1x,a,i5,2x,a3,1x,a,i5)')
     &                          nhelix,nhelix,rnam1,blah(ijk:ijk),ibeg,
     &                                        rnam2,blah(ijk:ijk),iend
                  i_in_helix=0
                endif
              endif
            enddo

C now do the sheets

            nsheet=0
            i_in_sheet=0
            do i=1,ires_max
              if(dssp_tmp(i:i).eq.'E') then
                if(i_in_sheet.eq.0) then
                  nsheet=nsheet+1
                  i_in_sheet=1
                  ibeg=i
                  iend=i
                elseif(i_in_sheet.eq.1) then
                  iend=i
                endif
              else
                if(i_in_sheet.eq.1) then
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'A') rnam1='ALA'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'C') rnam1='CYS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'D') rnam1='ASP'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'E') rnam1='GLU'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'F') rnam1='PHE'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'G') rnam1='GLY'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'H') rnam1='HIS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'I') rnam1='ILE'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'K') rnam1='LYS'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'L') rnam1='LEU'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'M') rnam1='MET'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'N') rnam1='ASN'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'P') rnam1='PRO'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'Q') rnam1='GLN'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'R') rnam1='ARG'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'S') rnam1='SER'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'T') rnam1='THR'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'V') rnam1='VAL'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'W') rnam1='TRP'
                if(seq_of_fasta(iii)(ibeg:ibeg).eq.'Y') rnam1='TYR'
                if(seq_of_fasta(iii)(iend:iend).eq.'A') rnam2='ALA'
                if(seq_of_fasta(iii)(iend:iend).eq.'C') rnam2='CYS'
                if(seq_of_fasta(iii)(iend:iend).eq.'D') rnam2='ASP'
                if(seq_of_fasta(iii)(iend:iend).eq.'E') rnam2='GLU'
                if(seq_of_fasta(iii)(iend:iend).eq.'F') rnam2='PHE'
                if(seq_of_fasta(iii)(iend:iend).eq.'G') rnam2='GLY'
                if(seq_of_fasta(iii)(iend:iend).eq.'H') rnam2='HIS'
                if(seq_of_fasta(iii)(iend:iend).eq.'I') rnam2='ILE'
                if(seq_of_fasta(iii)(iend:iend).eq.'K') rnam2='LYS'
                if(seq_of_fasta(iii)(iend:iend).eq.'L') rnam2='LEU'
                if(seq_of_fasta(iii)(iend:iend).eq.'M') rnam2='MET'
                if(seq_of_fasta(iii)(iend:iend).eq.'N') rnam2='ASN'
                if(seq_of_fasta(iii)(iend:iend).eq.'P') rnam2='PRO'
                if(seq_of_fasta(iii)(iend:iend).eq.'Q') rnam2='GLN'
                if(seq_of_fasta(iii)(iend:iend).eq.'R') rnam2='ARG'
                if(seq_of_fasta(iii)(iend:iend).eq.'S') rnam2='SER'
                if(seq_of_fasta(iii)(iend:iend).eq.'T') rnam2='THR'
                if(seq_of_fasta(iii)(iend:iend).eq.'V') rnam2='VAL'
                if(seq_of_fasta(iii)(iend:iend).eq.'W') rnam2='TRP'
                if(seq_of_fasta(iii)(iend:iend).eq.'Y') rnam2='TYR'
                write(32,'("SHEET",i5,3x,"A",1x,"1",1x,a3,1x,a,i4,2x,
     &                                                 a3,1x,a,i4)')
     &                                 nsheet,rnam1,blah(ijk:ijk),ibeg,
     &                                        rnam2,blah(ijk:ijk),iend
                write(22,'("SHEET",i5,3x,"A",1x,"1",1x,a3,1x,a,i4,2x,
     &                                                 a3,1x,a,i4)')
     &                                 nsheet,rnam1,blah(ijk:ijk),ibeg,
     &                                        rnam2,blah(ijk:ijk),iend
                  i_in_sheet=0
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      close(93)


C write out a reminder of the beta factor definitions

      write(32,'("REMARK:")')
      write(32,'("REMARK:   bfactor =  0.00 : model-built residue")')
      write(32,'("REMARK:   bfactor = 50.00 : aligned but dissimilar")')
      write(32,'("REMARK:   bfactor = 75.00 : aligned and similar")')
      write(32,'("REMARK:   bfactor = 99.99 : aligned and identical")')
      write(32,'("REMARK:")')
      write(32,'("REMARK:   occpncy =  0.00 : model-built residue")')
      write(32,'("REMARK:   occpncy =  1.00 : structured domain #1")')
      if(my_domain_max.ge.2) 
     &  write(32,'("REMARK:   occpncy =  2.00 : structured domain #2")')
      if(my_domain_max.ge.3) 
     &  write(32,'("REMARK:   occpncy =  3.00 : structured domain #3")')
      if(my_domain_max.ge.4) 
     &  write(32,'("REMARK:   occpncy =  4.00 : structured domain #4")')
      if(my_domain_max.ge.5) 
     &  write(32,'("REMARK:   occpncy =  5.00 : structured domain #5")')
      write(32,'("REMARK:")')

C come here to just write out the atoms of this model

4321  continue
      write(32,'("MODEL ",i8)')nmodels
      write(22,'("MODEL ",i8)')nmodels

C now start reading the last loopy pdb file for all coordinates

711   read(11,21,end=714)strong

C we need to skip any HETATM entries that were added at the loopy stage
C to prevent loops going crazy places - these entries can be
C distinguished by a -99.99 in the occupancy field - they will, of
C course, be added back later anyway...

      if(strong(55:60).eq.'-99.99') goto 711

C v4.5 we may have HETATM entries that don't conform to the above so
C let's also explicitly check for HETATM...

      if(strong(1:4).eq.'HETA') goto 711

C note that when we write out atoms we ask not only was it in the
C original pdb file (true for any kind of match), we also ask if it was
C model-built previously (i.e. it has a beta factor of zero) - if the
C latter is true we make sure it still gets beta=0 even though
C i_am_original=1 for it...

C UPDATE - we haven't yet figured out all the options here :(

      nats=nats+1
      read(strong,'(22x,i4)')ires
      read(strong,'(21x,a)')chainID
      read(strong,'(30x,3f8.3)')xfin,yfin,zfin ! use last cords for DUM

C store the coords so we can omit all "O   DUM/N   DUM" atoms from OPM
C that are too close to the protein - hopefully this will help avoid VMD
C showing truncated CARTOON representations...

      xlast(nats)=xfin
      ylast(nats)=yfin
      zlast(nats)=zfin

C v2.1 we have ires but we need to know the final chainID and the final
C res# - the easiest way to do this is loop over all fragments till we
C find the one that has the same value of cur_chainID and for which the
C first and last residue numbers encompass the just-read residue #...

C for cases where there is only one fragment set defaults:

c     if(num_frags.eq.1) then
c       ires_final=ires
c       chainID_final=fin_chainID(1)
c     else

C v2.9 we also find ires_total - i.e. res# in total-res system...

      ires_total=0
      ifound=0
      do nyyy=1,num_frags
        if(ires.gt.cur_offset(nyyy).and.
     &     ires.le.cur_offset(nyyy)+num_res_in_frag(nyyy).and.
     &     chainID.eq.cur_chainID(nyyy)) then
          my_frag=nyyy
          ires_final=ires-cur_offset(nyyy)
          chainID_final=fin_chainID(nyyy)
          ifound=1
          ires_total=ires_total+ires_final ! v2.9
          exit
        else
          ires_total=ires_total+num_res_in_frag(nyyy) ! v2.9
        endif
      enddo
      if(ifound.eq.0) then
        write(*,*)
        write(*,*)'shit we did not find ires = ',ires
        write(*,*)'and chainID = ',chainID
        write(*,*)'when about to write out - this is fatal :('
        write(*,*)
        stop
      endif
c     endif

C write out the domain # to the occupancy field
C and the i_am_original setting to the beta factor field

      write(32,812)strong(1:4),nats,strong(12:21),
Cv2.1&             chn_tot_res(ires:ires),
Cv2.1&             res_tot_res(ires),strong(27:54),
     &             chainID_final,ires_final,strong(27:54), ! v2.1

C v2.9 - use next line if you want domain# written to occupancy
C v3.2 - bugfix here and below - use ires_total, not ires

Cv2.9&             real(my_domain(ires)),my_bfactr(ires)
     &             real(my_domain(ires_total)),my_bfactr(ires_total),

C v2.9 - use next line if you want psipred value written to occupancy

Cv2.9&             real(occ_psipred(ires_total)),my_bfactr(ires)
C    &             real(occ_psipred(ires_total)),my_bfactr(ires_total),

C v3.2 add on the total chain # for easy sanity checking

     &             my_frag

812   format(a4,i7,a10,a,i4,a28,2f6.2,' tot_chain#= ',i6)

C also write out to CARTOON if appropriate

      if(strong(13:16).eq.' N  '.or.
     &   strong(13:16).eq.' CA '.or.
     &   strong(13:16).eq.' C  '.or.
     &   strong(13:16).eq.' O  ') then
        nats2=nats2+1
        write(22,812)strong(1:4),nats2,strong(12:21),
Cv2.1&               chn_tot_res(ires:ires),
Cv2.1&               res_tot_res(ires),strong(27:54),
     &               chainID_final,ires_final,strong(27:54), ! v2.1

C v2.9 - use next line if you want domain# written to occupancy
C v3.2 - bugfix here and below - use ires_total, not ires

Cv2.9&               real(my_domain(ires)),my_bfactr(ires)
C    &               real(my_domain(ires_total)),my_bfactr(ires_total),

C v2.9 - use next line if you want psipred value written to occupancy

Cv2.9&               real(occ_psipred(ires_total)),my_bfactr(ires)
     &             real(occ_psipred(ires_total)),my_bfactr(ires_total),

C v3.2 add on the total chain # for easy sanity checking

     &               my_frag

      endif

C following lines are to write out structured values for residues that
C were not built by loopy here and that were not previously unstructd

c     elseif(inot_made_yet(ires).eq.0.and.
c    &       i_am_by_loopy(ires).eq.0.and. 
c    &       i_old_unstruc(ires).eq.0) then
c       write(32,712)strong(1:4),nats,strong(12:21),
c    &               chn_tot_res(ires:ires),
c    &               res_tot_res(ires),strong(27:60)

ctt   elseif(i_done_res(j).eq.1) then
ctt     write(32,712)strong(1:4),nats,strong(12:21),
ctt  &               chn_tot_res(ires:ires),
ctt  &               res_tot_res(ires),strong(27:60)

C AHE notes - I'm not sure that the following would ever be executed -
C not sure how something could match both criteria...

c     elseif(i_am_original(ires).ge.1.and.
c    &       i_old_unstruc(ires).eq.1) then
c       write(32,713)strong(1:4),nats,strong(12:21),
c    &               chn_tot_res(ires:ires),
c    &               res_tot_res(ires),strong(27:60)
c
c     else
c       write(32,713)strong(1:4),nats,strong(12:21),
c    &               chn_tot_res(ires:ires),
c    &               res_tot_res(ires),strong(27:60)
713     format(a4,i7,a10,a,i4,a34,'  0.00')
c     endif

C increment total charge on structure

      if(strong(13:16).eq.' N  '.and.
     &   strong(18:20).eq.'ARG') then
        qqq_tot_ahemodel=qqq_tot_ahemodel+1.0
      elseif(strong(13:16).eq.' N  '.and.
     &   strong(18:20).eq.'LYS') then
        qqq_tot_ahemodel=qqq_tot_ahemodel+1.0
      elseif(strong(13:16).eq.' N  '.and.
     &   strong(18:20).eq.'HIS') then
        qqq_tot_ahemodel=qqq_tot_ahemodel+0.5
      elseif(strong(13:16).eq.' N  '.and.
     &   strong(18:20).eq.'ASP') then
        qqq_tot_ahemodel=qqq_tot_ahemodel-1.0
      elseif(strong(13:16).eq.' N  '.and.
     &   strong(18:20).eq.'GLU') then
        qqq_tot_ahemodel=qqq_tot_ahemodel-1.0
      endif

      goto 711
714   close(11)

C before finally finally finishing we should add on the HETERO atoms
C note that we set all b-factors to zero so that they don't mess up the
C scale in VMD... we write these atoms to both AHEMODEL & CARTOON.pdbs

      nats_original=nats
      nhats=0 ! v4.5
      open(unit=21,file='HETERO_original.pdb',status='unknown')
171   read(21,21,end=173)strong
      strong(61:66)='  0.00'

C v4.5 increment # of HETATMs found and then ask if this is HETATM that
C got taken care of by a tail - if so, it'll have mh=0 and we skip it
C UPDATE - use simple 0/1 flag ih() instead...

      nhats=nhats+1
      if(ih(nhats).eq.0) goto 171

C v4.6 separate paths depending on whether we want to separate HETATMs
C if not, then just write out HETATMs to AHEMODEL and CARTOON pdbs
C if separating, then only write to HETATMs.pdb on the first nmodel...

      if(.not.separate_hetatms) then ! default prior to v4.6

        nats=nats+1
        write(32,715)strong(1:6),mod(nats,100000),strong(12:80)
        nats2=nats2+1
        write(22,715)strong(1:6),mod(nats2,100000),strong(12:80)

      elseif(separate_hetatms) then

        if(nmodels.eq.1) then
          nats=nats+1
          write(23,715)strong(1:6),mod(nats,100000),strong(12:80)
        endif

      endif

715   format(a6,i5,a69)
      goto 171
173   close(21)

C v4.6 if separate_hetatms then we can close the file...

      if(separate_hetatms) then
        if(nmodels.eq.1) then
          write(*,*)
          write(*,*)'finished writing fixed HETATMs, now closing '
          write(*,*)'file = ',sys_string9c(1:mmm3)
          write(*,*)
          close(23)
        endif
      endif

C v4.5 we should add on tail-associated HETATM entries here - these
C should have been skipped above - if not then we'll erroneously have
C two copies of the same HETATM entry: the original version, which will
C have weird coordinates, and the new version which is tail-attached...

C v4.6 - we need to make sure that these are out in chain order - in
C loopy_inp.pdb these HETATMs will be ordered according to the order in
C which the tails were built - so we first need to go through and assign
C a chain# (based on the chainID) and then write out in chain order...
C so the first time through we don't write anything out - we just
C figure out what chain#s we are having to deal with - note that this
C will fail if num_frags>78...

      ncrap=0
      ichain=0 
      open(unit=11,file='loopy_inp.pdb',status='unknown')
891   read(11,'(a66)',end=892)char66
      if(char66(1:4).ne.'HETA') goto 891

ccc   char66(61:66)='  0.00' ! again set b-factor to zero (see above)
      ncrap=ncrap+1
      do nyyy=1,num_frags
        if(char66(22:22).eq.blah(nyyy:nyyy)) then
          ichain(ncrap)=nyyy
          exit
        endif
      enddo
      if(ichain(ncrap).eq.0) then
        write(*,*)
        write(*,*)'failed to find the frag# for a HETATM in '
        write(*,*)'loopy_inp.pdb with a chainID ',char66(22:22)
        write(*,*)'this is a fatal error - may mean you have more'
        write(*,*)'than 78 chains...'
        write(*,*)
        stop
      endif

ccc   if(char66(1:4).eq.'HETA') write(32,'(a66)') char66
ccc   if(char66(1:4).eq.'HETA') write(22,'(a66)') char66
      goto 891
892   rewind(11)

C figure out max chain value and min chain value

      ichain_max=-1
      ichain_min=999999
      do nz=1,ncrap
        ichain_max=max(ichain_max,ichain(nz))
        ichain_min=min(ichain_min,ichain(nz))
      enddo

C now cycle from ichain_min to ichain_max reading loopy_inp.pdb again
C and again, each time writing out more HETATMs...

      do nz=ichain_min,ichain_max
        write(*,*)'finding domain-associated HETATMs for chain# ',nz
        ncrap=0
893     read(11,'(a66)',end=894)char66
        if(char66(1:4).ne.'HETA') goto 893
        ncrap=ncrap+1
        if(ichain(ncrap).ne.nz) goto 893
        char66(61:66)='  0.00' ! again set b-factor to zero (see above)
        write(32,'(a66)') char66
        write(22,'(a66)') char66
        goto 893
894     rewind(11)
      enddo
      close(11)

C now add on any OPM-DUM atoms that are not clashing with other atoms -
C note that we only write these out to CARTOON.pdb...
C v2.9 - only do this if membrane_dive!!!

      if(membrane_dive) then
        open(unit=21,file='DUM_original.pdb',status='unknown')
271     read(21,21,end=273)strong
        ikeep=1 ! assume we'll keep it till told otherwise
        read(strong,'(30x,3f8.3)')xcur,ycur,zcur
        do n=1,nats_original
          dist2=(xlast(n)-xcur)**2+
     &          (ylast(n)-ycur)**2+
     &          (zlast(n)-zcur)**2
          if(dist2.le.25.0) then
            ikeep=0 ! don't keep it if its <5A from a protein atom
            exit
          endif
        enddo
        if(ikeep.eq.0) then
          write(*,*)'skipping DUM # ',strong(13:26)
          goto 271 
        endif
        nats2=nats2+1
        strong(61:66)='  0.00'
        write(22,715)strong(1:6),mod(nats2,100000),strong(12:80)
        goto 271
273     close(21)
      endif

C finally, write out two dummy atoms with extreme beta values to set the
C scale properly in VMD - this covers all possible scenarios... note
C that we only write these out to CARTOON.pdb...
C v2.9 - we now also do this for occupancy so that colors look right for
C PSIPRED predictions - values of -1.0 and +1.0 are the extremes

      betastring( 1:30)='HETATM XXXX  DUM DUM     1   '
      nats2=nats2+1
      write(betastring(7:12),'(i5)')nats2
      write(betastring(31:54),'(3f8.3)')xfin,yfin,zfin
      betastring(55:66)=' -1.00  0.00'
      write(22,'(a66)')betastring

      nats2=nats2+1
      write(betastring(7:12),'(i5)')nats2
      write(betastring(31:54),'(3f8.3)')xfin,yfin,zfin
      betastring(55:66)='  1.00 99.99'
      write(22,'(a66)')betastring
      close(21)

      write(32,'("ENDMDL")')
      write(22,'("ENDMDL")')

C V3.7 - just a sanity check here: all models should have the same
C number of atoms...

      write(*,*)'FOR MODEL # ',nmodels,' WE WROTE ',nats,' ATOMS'

C v4.5 9 April 2021 reset ih using ih_original for all HETATM entries 
C we need to do this for each new MODEL that we build

        do n=1,num_het_atms
          ih(n)=ih_original(n)
        enddo

      if(nmodels.lt.nmodels_desired) goto 1234

C now that we have enough models we can finally close AHEMODEL and
C CARTOON pdb files...
     
      close(32)
      close(22)

      write(*,939)qqq_tot_scwrl,qqq_tot_ahemodel
939   format('Net Charge Before =  ',f15.3,' After = ',f15.5)
      write(*,941)
941   format('these charges are only approximate (+0.5 for His etc)')

C finally, clean up the temp4_files...

      call system('rm temp4_????_????')

C v4.7 and finally, finally, evaluate the total time taken:

      call system_clock(count=ktime2,count_rate=ru)
      time_total=real(ktime2-ktime1)/real(ru)

      write(*,*)
      write(*,*)'time_total                 = ',time_total
      write(*,*)

      stop
      end



      subroutine matrix(om1,om2,om3,rot11,rot12,rot13,
     &                              rot21,rot22,rot23,
     &                              rot31,rot32,rot33)
      implicit real(a-h,o-z)
      c1=cos(om1)
      s1=sin(om1)
      c2=cos(om2)
      s2=sin(om2)
      c3=cos(om3)
      s3=sin(om3)
      rot11=c2*c3
      rot12=s1*s2*c3-c1*s3
      rot13=s1*s3+c1*s2*c3
      rot21=c2*s3
      rot22=c1*c3+s1*s2*s3
      rot23=c1*s2*s3-s1*c3
      rot31=-s2
      rot32=s1*c2
      rot33=c1*c2
      return
      end


      real function ran1(idumahe)
C
C note that this is the new ran routine from numerical receipes for F90
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: K4B=selected_int_kind(9)
      INTEGER(K4B), INTENT(INOUT) :: idumahe
      REAL :: ran
      INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,
     &                           IQ=127773,IR=2836
      REAL, SAVE :: am
      INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
      if (idumahe <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/IM
        iy=ior(ieor(888889999,abs(idumahe)),1)
        ix=ieor(777755555,abs(idumahe))
        idumahe=abs(idumahe)+1
      endif
      ix=ieor(ix,ishft(ix,13))
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/IQ
      iy=IA*(iy-k*IQ)-IR*k
      if (iy < 0) iy=iy+IM
      ran1=am*ior(iand(IM,ieor(ix,iy)),1)
      END FUNCTION ran1



      subroutine acrosb (a,b,c)
      implicit real (a-h,o-z)
      dimension a(3),b(3),c(3)
C
C-----------------------------------------------------------------------
C
C      c=a x b - c returned normalized
C
C-----------------------------------------------------------------------
C
      x = a(2)*b(3)-b(2)*a(3)
      y = b(1)*a(3)-a(1)*b(3)
      z = a(1)*b(2)-b(1)*a(2)
      s = 1.000/sqrt(x**2+y**2+z**2)
      c(1) = x*s
      c(2) = y*s
      c(3) = z*s
      return
      end



      subroutine pdbsup(nat_loc,xm_loc,ym_loc,zm_loc,rf_loc,rms_loc)

C note that nat_loc is the actual number of atoms in this molecule+4
C with the 4 being the N,CA,C,O atoms...

c ----------------------------------------------------------------------PDBS0002
c                               PROGRAM PDBSUP                          PDBS0001
c             Determines rotation matrix and translation vector for     PDBS0002
c              best fit superimposition of two pdb files by solving     PDBS0003
c                      the quaternion eigenvalue problem.               PDBS0004
c                      Code by B.Rupp and S.Parkin (1996)               PDBS0005
c               Lawrence Livermore National Laboratory, BBRP            PDBS0006
c             Method by S.K.Kearsley, Acta Cryst. A45, 208 (1989)       PDBS0007
c              See http://www-structure.llnl.gov for details.           PDBS0008
c ----------------------------------------------------------------------PDBS0009
c !MS$DEBUG                                                             PDBS0010
                                                                        PDBS0011
c --- number of atoms read -                                            PDBS0012

      implicit real(a-h,o-z)
                                                                        PDBS0014
      integer, intent(in)     :: nat_loc
      real, intent(inout)     :: xm_loc(1:nat_loc)
      real, intent(inout)     :: ym_loc(1:nat_loc)
      real, intent(inout)     :: zm_loc(1:nat_loc)
      real, intent(in)        :: rf_loc(1:4,1:3)
      real,    intent(out)    :: rms_loc

C note that we make rm_loc from xm_loc,ym_loc,zm_loc...

      real                    :: rm_loc(1:nat_loc,1:3)

      real cack(3)
      real dm(4),vm(4,4),cm(3),cf(3)                                    PDBS0018
      real tr(3),t(3,3),q(4,4)                                          PDBS0019
      real dxp(nat_loc,3),dxm(nat_loc,3)
      real sumf(3),summ(3)                                               PDBS0021
                                                                        PDBS0022
      deg = 57.29577951                                                 PDBS0023
      zero=10E-30                                                       PDBS0024

C AHE adds this to avoid having to reallocate rm every time...

      do j=1,nat_loc
        rm_loc(j,1)=xm_loc(j)
        rm_loc(j,2)=ym_loc(j)
        rm_loc(j,3)=zm_loc(j)
      enddo
        
      imov=4 ! use only 4 atoms for the superpositions... 

c --- initialize all --                                                 PDBS0124
      do i=1,3                                                          PDBS0125
         sumf(i)=0                                                      PDBS0126
         summ(i)=0                                                      PDBS0127
      end do                                                            PDBS0128
      call filmat(4,4,q,0)                                              PDBS0129
                                                                        PDBS0130
c --- sum up all coordinates (in dble precision) to find centre ---     PDBS0131
      do k=1,imov                                                       PDBS0132
         do i=1,3                                                       PDBS0133
            sumf(i)=sumf(i)+rf_loc(k,i)                                     PDBS0134
            summ(i)=summ(i)+rm_loc(k,i)                                     PDBS0135
         end do                                                         PDBS0num_solute_typs
      end do                                                            PDBS0137
      do i=1,3                                                          PDBS0138
         cm(i)=summ(i)/real(imov)                                           PDBS0139
         cf(i)=sumf(i)/real(imov)                                           PDBS0140
         tr(i)=cf(i)-cm(i)                                              PDBS0141
      end do                                                            PDBS0142
                                                                        PDBS0143
c     write(*,'(/a,3f12.3)')' Cen of target molecule  =',(cf(i),i=1,3)  PDBS0144
c     write(*,'(a,3f12.3)') ' Cen of moving molecule  =',(cm(i),i=1,3)  PDBS0145
c     write(*,'(a,3f8.3/)')' T - vector probe -> target =',(tr(i),i=1,3)PDBS0146
                                                                        PDBS0147
c     write (*,'(a)')' Creating coordinate differences.......'          PDBS0148
c --- create coordinate differences delta x plus (dxp) and minus (dxm)  PDBS0149
      do k=1,imov                                                       PDBS0150
         do j=1,3                                                       PDBS0151
            dxm(k,j)=rm_loc(k,j)-cm(j)-(rf_loc(k,j)-cf(j))                      PDBS0152
            dxp(k,j)=rm_loc(k,j)-cm(j)+(rf_loc(k,j)-cf(j))                      PDBS0153
         end do                                                         PDBS0154
      end do                                                            PDBS0155
                                                                        PDBS0156
c --- fill upper triangle of (symmetric) quaternion matrix --           PDBS0157
c     write (*,'(a)')' Filling quaternion matrix ............'          PDBS0158
      do k=1,imov                                                       PDBS0159
c ---    diags are sums of squared cyclic coordinate differences        PDBS0160
         q(1,1)=q(1,1)+dxm(k,1)**2+dxm(k,2)**2+dxm(k,3)**2              PDBS0161
         q(2,2)=q(2,2)+dxp(k,2)**2+dxp(k,3)**2+dxm(k,1)**2              PDBS0162
         q(3,3)=q(3,3)+dxp(k,1)**2+dxp(k,3)**2+dxm(k,2)**2              PDBS0163
         q(4,4)=q(4,4)+dxp(k,1)**2+dxp(k,2)**2+dxm(k,3)**2              PDBS0164
c ---    cross differences                                              PDBS0165
         q(1,2)=q(1,2)+dxp(k,2)*dxm(k,3)-dxm(k,2)*dxp(k,3)              PDBS0166
         q(1,3)=q(1,3)+dxm(k,1)*dxp(k,3)-dxp(k,1)*dxm(k,3)              PDBS0167
         q(1,4)=q(1,4)+dxp(k,1)*dxm(k,2)-dxm(k,1)*dxp(k,2)              PDBS0168
         q(2,3)=q(2,3)+dxm(k,1)*dxm(k,2)-dxp(k,1)*dxp(k,2)              PDBS0169
         q(2,4)=q(2,4)+dxm(k,1)*dxm(k,3)-dxp(k,1)*dxp(k,3)              PDBS0170
         q(3,4)=q(3,4)+dxm(k,2)*dxm(k,3)-dxp(k,2)*dxp(k,3)              PDBS0171
      end do                                                            PDBS0172
c --- fill the rest by transposing it onto itself                       PDBS0173
      call trpmat(4,q,q)                                                PDBS0174
c     write (*,'(/a)')                                                  PDBS0175
c    &     '       q(1)         q(2)         q(3)        q(4)'          PDBS0176
c     do i=1,4                                                          PDBS0177
c        write(*,'(4e13.5)') (q(i,j),j=1,4)                             PDBS0178
c     end do                                                            PDBS0179
                                                                        PDBS0180
c --- orthogonalization by jacobi rotation = solution of EV -problem -- PDBS0181
c     write (*,'(/a)')' Jacobi orthogonalization ..........'            PDBS0182
      n=4                                                               PDBS0183
      ns=4                                                              PDBS0184
      call jacobi(q,n,ns,dm,vm,nmrot)                                   PDBS0185
c --- sort eigenvectors after eigenvalues, descending --                PDBS0186
c     write (*,'(a/)')' Sorting eigenvalues/vectors .......'            PDBS0187
      call eigsrt(dm,vm,n,ns)                                           PDBS0188
c     write (*,'(a,i2,a)')' Eigenvalues and Eigenvectors (',            PDBS0189
c    & nmrot,' Jacobi rotations)'                                       PDBS0190
c     write (*,'(a)') '      e(1)        e(2)        e(4)        e(4)'  PDBS0191
c     write (*,'(4e12.5,i5)') (dm(j),j=1,4)                             PDBS0192
c     write (*,'(a)') '      ev(1)       ev(2)       ev(3)       ev(4)' PDBS0193
c     do i=1,4                                                          PDBS0194
c        write(*,'(4f12.6)') (vm(i,j),j=1,4)                            PDBS0195
c     end do                                                            PDBS0196
                                                                        PDBS0197
c --- the smallest eigenvector contains best fit srs                    PDBS0198
      rmsd=sqrt(abs(dm(4)/imov))                                        PDBS0199

c     write(*,*)'my rmsd = ',rmsd
c     write(*,'(/a/)')                                                  PDBS0200
c    & ' The smallest eigenvalue represents s.r.s. of best fit'         PDBS0201
c     write(*,'(a)')                                                    PDBS0202
c    & ' Constructing the best fit rotation matrix from associated'     PDBS0203
c     write(*,'(a/)') ' eigenvector elements (last column).....'        PDBS0204
                                                                        PDBS0205
c --- fill the rotation matrix which is made of elements from 4th EV    PDBS0206
      t(1,1)=vm(1,4)**2+vm(2,4)**2-vm(3,4)**2-vm(4,4)**2                PDBS0207
      t(2,1)=2*(vm(2,4)*vm(3,4)+vm(1,4)*vm(4,4))                        PDBS0208
      t(3,1)=2*(vm(2,4)*vm(4,4)-vm(1,4)*vm(3,4))                        PDBS0209
      t(1,2)=2*(vm(2,4)*vm(3,4)-vm(1,4)*vm(4,4))                        PDBS0210
      t(2,2)=vm(1,4)**2+vm(3,4)**2-vm(2,4)**2-vm(4,4)**2                PDBS0211
      t(3,2)=2*(vm(3,4)*vm(4,4)+vm(1,4)*vm(2,4))                        PDBS0212
      t(1,3)=2*(vm(2,4)*vm(4,4)+vm(1,4)*vm(3,4))                        PDBS0213
      t(2,3)=2*(vm(3,4)*vm(4,4)-vm(1,4)*vm(2,4))                        PDBS0214
      t(3,3)=vm(1,4)**2+vm(4,4)**2-vm(2,4)**2-vm(3,4)**2                PDBS0215
                                                                        PDBS0216
c     do i=1,3                                                          PDBS0217
c        write(*,'(3f11.5)') (t(i,j),j=1,3)                             PDBS0218
c     end do                                                            PDBS0219
                                                                        PDBS0220
c --- reset dxm to store the individual rmsd's in it now -              PDBS0221
      call filmat(nat_loc,3,dxm,0)                                          PDBS0222
                                                                        PDBS0223
c --- xm and xf are not translated                                      PDBS0224
      do k=1,nat_loc                                                    PDBS0225
c ---    subtract cm                                                    PDBS0226
         rm_loc(k,1)=rm_loc(k,1)-cm(1)                                             PDBS0227
         rm_loc(k,2)=rm_loc(k,2)-cm(2)                                             PDBS0227
         rm_loc(k,3)=rm_loc(k,3)-cm(3)                                             PDBS0227
C
C AHE
C make call to rotvec easier
C
         cack(1)=rm_loc(k,1)
         cack(2)=rm_loc(k,2)
         cack(3)=rm_loc(k,3)
c ---    rotate it                                                      PDBS0228
         call rotvec(3,cack,t)                                          PDBS0229
c ---    now add cf                                                     PDBS0230
         xm_loc(k)=cack(1)
         ym_loc(k)=cack(2)
         zm_loc(k)=cack(3)
         xm_loc(k)=xm_loc(k)+cf(1)                                          PDBS0231
         ym_loc(k)=ym_loc(k)+cf(2)                                          PDBS0231
         zm_loc(k)=zm_loc(k)+cf(3)                                          PDBS0231
c        do i=1,3                                                       PDBS0232
c           dxm(k,i)=sqrt((rf_loc(k,i)-rm_loc(k,i))**2)                         PDBS0233
c        end do                                                         PDBS0234
      end do                                                            PDBS0235
                                                                        PDBS0236
      rms_loc=rmsd
c     write(6,3)s                                                       PDBS0260
                                                                        PDBS0293
 0003      format(a)

      end                                                               PDBS0354
                                                                        PDBS0355
      subroutine inkey (answ)                                           INKE0001
c ----------------------------------------------------------------------INKE0002
c     reads a key as answer                                             INKE0003
c ----------------------------------------------------------------------INKE0004
      character answ                                                    INKE0005
      read (*,'(a1)') answ                                              INKE0006
      call upstrg (answ,1)                                              INKE0007
      if ((answ.eq.' ').or.(answ.eq.char(13))) answ='Y'                 INKE0008
      return                                                            INKE0009
      end                                                               INKE0010
                                                                        INKE0011
      subroutine upstrg(strg,istrln)                                    UPST0001
c ----------------------------------------------------------------------UPST0002
c converts string str$ of lenght istrlen to upcase                      UPST0003
c ----------------------------------------------------------------------UPST0004
      character strg(istrln)                                            UPST0005
      integer      iascii                                               UPST0006
C-----change to upper case                                              UPST0007
         do 3031 i=1,istrln                                             UPST0008
            if (ichar(strg(i)).ge.95.and.ichar(strg(i)).le.122) then    UPST0009
               iascii=ichar(strg(i))-32                                 UPST0010
               strg(i)=char(iascii)                                     UPST0011
            endif                                                       UPST0012
 3031    continue                                                       UPST0013
      return                                                            UPST0014
      end                                                               UPST0015
                                                                        UPST0016
      subroutine trpmat(n,t,tr)                                         TRPM0001
c --- transpose matrix -------------------------------------------------TRPM0002
      real t(n,n), tr(n,n)                                              TRPM0003
      do i=1,n                                                          TRPM0004
         do j=1,n                                                       TRPM0005
            tr(j,i)=t(i,j)                                              TRPM0006
         end do                                                         TRPM0007
      end do                                                            TRPM0008
      return                                                            TRPM0009
      end                                                               TRPM0010
                                                                        TRPM0011
      subroutine filmat(n,m,r,ifil)                                     FILM0001
      real r(n,m)                                                       FILM0002
c --- initialize matrix ------------------------------------------------FILM0003
      do i=1,n                                                          FILM0004
         do j=1,m                                                       FILM0005
            r(i,j)=ifil                                                 FILM0006
         end do                                                         FILM0007
      end do                                                            FILM0008
      return                                                            FILM0009
      end                                                               FILM0010
                                                                        FILM0011
      subroutine rotvec (n,v,t)                                         ROTV0001
c --- multiply vector with matrix --------------------------------------ROTV0002
      real t(n,n), v(n),s(n)                                            ROTV0003
                                                                        ROTV0004
      do i=1,n                                                          ROTV0005
         s(i)=v(i)                                                      ROTV0006
         v(i)=0.0                                                       ROTV0007
      end do                                                            ROTV0008
      do i=1,n                                                          ROTV0009
         do j=1,n                                                       ROTV0010
            v(i)=v(i)+s(j)*t(i,j)                                       ROTV0011
         end do                                                         ROTV0012
      end do                                                            ROTV0013
      return                                                            ROTV0014
      end                                                               ROTV0015
                                                                        ROTV0016
      SUBROUTINE eigsrt(d,v,n,np)                                       EIGS0001
c ----------------------------------------------------------------------EIGS0002
      INTEGER n,np                                                      EIGS0003
      REAL d(np),v(np,np)                                               EIGS0004
      INTEGER i,j,k                                                     EIGS0005
      REAL p                                                            EIGS0006
      do 13 i=1,n-1                                                     EIGS0007
        k=i                                                             EIGS0008
        p=d(i)                                                          EIGS0009
        do 11 j=i+1,n                                                   EIGS0010
          if(d(j).ge.p)then                                             EIGS0011
            k=j                                                         EIGS0012
            p=d(j)                                                      EIGS0013
          endif                                                         EIGS0014
11      continue                                                        EIGS0015
        if(k.ne.i)then                                                  EIGS0016
          d(k)=d(i)                                                     EIGS0017
          d(i)=p                                                        EIGS0018
          do 12 j=1,n                                                   EIGS0019
            p=v(j,i)                                                    EIGS0020
            v(j,i)=v(j,k)                                               EIGS0021
            v(j,k)=p                                                    EIGS0022
12        continue                                                      EIGS0023
        endif                                                           EIGS0024
13    continue                                                          EIGS0025
      return                                                            EIGS0026
      END                                                               EIGS0027
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             EIGS0028
                                                                        EIGS0029
      SUBROUTINE jacobi(a,n,np,d,v,nrot)                                JACO0001
c ----------------------------------------------------------------------JACO0002
c     modified from numerical recipes book                              JACO0003
c     one needs to set the threshold for sm from sm.eq.0 to sm.lt.10E-30JACO0004
c     (anything in this range would be ok) due to underflow errors on   JACO0005
c     some computers/compilers.                                         JACO0006
c ----------------------------------------------------------------------JACO0007
      PARAMETER (nmax=500)                                              JACO0008
                                                                        JACO0009
      INTEGER n,np,nrot                                                 JACO0010
      REAL a(np,np),d(np),v(np,np)                                      JACO0011
      INTEGER i,ip,iq,j,maxrot                                          JACO0012
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),zero          JACO0013
                                                                        JACO0014
c --- zero set and iteration maximum                                    JACO0015
      zero=10E-30                                                       JACO0016
      maxrot=50                                                         JACO0017
                                                                        JACO0018
      do 12 ip=1,n                                                      JACO0019
        do 11 iq=1,n                                                    JACO0020
          v(ip,iq)=0.                                                   JACO0021
11      continue                                                        JACO0022
        v(ip,ip)=1.                                                     JACO0023
12    continue                                                          JACO0024
      do 13 ip=1,n                                                      JACO0025
        b(ip)=a(ip,ip)                                                  JACO0026
        d(ip)=b(ip)                                                     JACO0027
        z(ip)=0.                                                        JACO0028
13    continue                                                          JACO0029
      nrot=0                                                            JACO0030
      do 24 i=1,maxrot                                                  JACO0031
        sm=0.                                                           JACO0032
        do 15 ip=1,n-1                                                  JACO0033
          do 14 iq=ip+1,n                                               JACO0034
            sm=sm+abs(a(ip,iq))                                         JACO0035
14        continue                                                      JACO0036
15      continue                                                        JACO0037
c ---   modified convergence threshold ---                              JACO0038
        if(sm.lt.zero)return                                            JACO0039
        if(i.lt.4)then                                                  JACO0040
          tresh=0.2*sm/n**2                                             JACO0041
        else                                                            JACO0042
          tresh=0.                                                      JACO0043
        endif                                                           JACO0044
        do 22 ip=1,n-1                                                  JACO0045
          do 21 iq=ip+1,n                                               JACO0046
            g=100.*abs(a(ip,iq))                                        JACO0047
            if((i.gt.4).and.(abs(d(ip))+                                JACO0048
     &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then            JACO0049
              a(ip,iq)=0.                                               JACO0050
            else if(abs(a(ip,iq)).gt.tresh)then                         JACO0051
              h=d(iq)-d(ip)                                             JACO0052
              if(abs(h)+g.eq.abs(h))then                                JACO0053
                t=a(ip,iq)/h                                            JACO0054
              else                                                      JACO0055
                theta=0.5*h/a(ip,iq)                                    JACO0056
                t=1./(abs(theta)+sqrt(1.+theta**2))                     JACO0057
                if(theta.lt.0.)t=-t                                     JACO0058
              endif                                                     JACO0059
              c=1./sqrt(1+t**2)                                         JACO0060
              s=t*c                                                     JACO0061
              tau=s/(1.+c)                                              JACO0062
              h=t*a(ip,iq)                                              JACO0063
              z(ip)=z(ip)-h                                             JACO0064
              z(iq)=z(iq)+h                                             JACO0065
              d(ip)=d(ip)-h                                             JACO0066
              d(iq)=d(iq)+h                                             JACO0067
              a(ip,iq)=0.                                               JACO0068
              do 16 j=1,ip-1                                            JACO0069
                g=a(j,ip)                                               JACO0070
                h=a(j,iq)                                               JACO0071
                a(j,ip)=g-s*(h+g*tau)                                   JACO0072
                a(j,iq)=h+s*(g-h*tau)                                   JACO0073
16            continue                                                  JACO0074
              do 17 j=ip+1,iq-1                                         JACO0075
                g=a(ip,j)                                               JACO0076
                h=a(j,iq)                                               JACO0077
                a(ip,j)=g-s*(h+g*tau)                                   JACO0078
                a(j,iq)=h+s*(g-h*tau)                                   JACO0079
17            continue                                                  JACO0080
              do 18 j=iq+1,n                                            JACO0081
                g=a(ip,j)                                               JACO0082
                h=a(iq,j)                                               JACO0083
                a(ip,j)=g-s*(h+g*tau)                                   JACO0084
                a(iq,j)=h+s*(g-h*tau)                                   JACO0085
18            continue                                                  JACO0086
              do 19 j=1,n                                               JACO0087
                g=v(j,ip)                                               JACO0088
                h=v(j,iq)                                               JACO0089
                v(j,ip)=g-s*(h+g*tau)                                   JACO0090
                v(j,iq)=h+s*(g-h*tau)                                   JACO0091
19            continue                                                  JACO0092
              nrot=nrot+1                                               JACO0093
            endif                                                       JACO0094
21        continue                                                      JACO0095
22      continue                                                        JACO0096
        do 23 ip=1,n                                                    JACO0097
          b(ip)=b(ip)+z(ip)                                             JACO0098
          d(ip)=b(ip)                                                   JACO0099
          z(ip)=0.                                                      JACO0100
23      continue                                                        JACO0101
24    continue                                                          JACO0102
ccc   pause 'too many iterations in jacobi'                             JACO0103
      return                                                            JACO0104
      END                                                               JACO0105
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             JACO0106



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


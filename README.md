## IPyRSSA

IPyRSSA (<b>I</b>ntegrative <b>Py</b>thon library for <b>R</b>NA <b>S</b>econdary <b>S</b>tructure <b>A</b>nalysis) is a set of Python library to analyze RNA secondary structure and SHAPE data.

#### New
Python 3 is supported now. Python 2 is not supported.

#### Update your local library

`git pull origin `

<h3> General module </h3>
`import General`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> load_fasta </td>
	<td> Read fasta file </td>
</tr>
<tr>
	<td> write_fasta </td>
	<td> Write fasta file </td>
</tr>
<tr>
	<td> load_dot </td>
	<td> Read dotBracket file </td>
</tr>
<tr>
	<td> write_dot </td>
	<td> Write dotBracket file </td>
</tr>
<tr>
	<td> load_shape </td>
	<td> Read SHAPE .out file </td>
</tr>
<tr>
	<td> load_SHAPEMap </td>
	<td> Read SHAPEmap file </td>
</tr>
<tr>
	<td> load_ct </td>
	<td> Read .ct file </td>
</tr>
<tr>
	<td> write_ct </td>
	<td> Write .ct file </td>
</tr>
<tr>
	<td> init_pd_rect </td>
	<td> Build a dataframe </td>
</tr>
<tr>
	<td> init_list_rect </td>
	<td> Build a list matrix </td>
</tr>
<tr>
	<td> find_all_match </td>
	<td> Find all match regions with a regex </td>
</tr>
<tr>
	<td> bi_search </td>
	<td> Binary search </td>
</tr>
<tr>
	<td> calc_shape_gini </td>
	<td> Calculate SHAPE gini index </td>
</tr>
<tr>
	<td> calc_shape_structure_ROC </td>
	<td> Calculate the ROC points structure and shape scores </td>
</tr>
<tr>
	<td> calc_AUC </td>
	<td> Calculate AUC with ROC points </td>
</tr>
<tr>
	<td> calc_AUC_v2 </td>
	<td> Calculate AUC with dot and shape_list </td>
</tr>
<tr>
	<td> seq_entropy </td>
	<td> Calculate the entropy of the sequence. </td>
</tr>
</table>


<h3> Colors module </h3>
`import Colors`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> format or f </td>
	<td> Format a colorful text </td>
</tr>
<tr>
	<td> color_SHAPE </td>
	<td> Convert SHAPE list to colorful blocks </td>
</tr>
<tr>
	<td> color_Seq_SHAPE </td>
	<td> Convert sequence to colorful sequence </td>
</tr>
<tr>
	<td> browse_shape </td>
	<td> Print and compare single/multiple shape scores <a href="examples/visual_icSHAPE_in_terminal">example</a> </td>
</tr>
<tr>
	<td> browse_multi_shape </td>
	<td> Align multiple sequences and print shape scores <a href="examples/visual_icSHAPE_in_terminal">example</a> </td>
</tr>
</table>

<h3> Cluster module </h3>
`import Cluster`

Warning: This module can only be used on loginviewxx/mgtxx

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> new_job </td>
	<td> Get a job handle </td>
</tr>
<tr>
	<td> handle.set_job_depends </td>
	<td> The job will be executed when parameter jobs done </td>
</tr>
<tr>
	<td> handle.submit </td>
	<td> Submit the job to queue </td>
</tr>
<tr>
	<td> handle.has_finish </td>
	<td> Return True if finished </td>
</tr>
<tr>
	<td> handle.job_status </td>
	<td> Return one of Not_Found, DONE, RUN, PEND, EXIT </td>
</tr>
<tr>
	<td> handle.wait </td>
	<td> Wait the job to finish </td>
</tr>
<tr>
	<td> handle.kill </td>
	<td> Kill the job </td>
</tr>
</table>

<h3> Seq module </h3>
`import Seq`

Prerequisites: pyliftover, pysam

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> reverse_comp </td>
	<td> Get reversed complementary sequence of raw sequence </td>
</tr>
<tr>
	<td> flat_seq </td>
	<td> Flatten the long sequence to multiline sequence </td>
</tr>
<tr>
	<td> format_gene_type </td>
	<td> Classify the raw gene type in annotation to common gene type </td>
</tr>
<tr>
	<td> Class:seqClass </td>
	<td> A class to fetch sequence from big genome </td>
</tr>
<tr>
	<td> lift_genome </td>
	<td> Convert the genome version (hg19=>hg38) </td>
</tr>
<tr>
	<td> search_subseq_from_genome </td>
	<td> Search a pattern in genome region </td>
</tr>
</table>


<h3> Structure module </h3>
`import Structure`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> predict_structure </td>
	<td> Prediction secondary structure combine SHAPE or not </td>
</tr>
<tr>
	<td> bi_fold </td>
	<td> Prediction RNA interaction </td>
</tr>
<tr>
	<td> search_TT_cross_linking </td>
	<td> Search TT cross linking sites in structure </td>
</tr>
<tr>
	<td> dyalign </td>
	<td> Predict a common secondary structure for two sequences </td>
</tr>
<tr>
	<td> multialign </td>
	<td> Predict a common secondary structure for multiple sequences </td>
</tr>
<tr>
	<td> estimate_energy </td>
	<td> Calculate the folding free energy change of a structure </td>
</tr>
<tr>
	<td> partition </td>
	<td> Calculate the partition function </td>
</tr>
<tr>
	<td> maxExpect </td>
	<td> Calculate the max-expect structure </td>
</tr>
<tr>
	<td> evaluate_dot </td>
	<td> Evaluate the Sensitivty and PPV for a predicted structure relative to target structure </td>
</tr>
<tr>
	<td> calc_structure_similarity </td>
	<td> Calculate the structure similarity,distance </td>
</tr>
<tr>
	<td> dot2ct </td>
	<td> Dotbracket to list </td>
</tr>
<tr>
	<td> dot2bpmap </td>
	<td> Dotbracket to dictionary </td>
</tr>
<tr>
	<td> parse_pseudoknot </td>
	<td> Parse pseudoknot with ctList </td>
</tr>
<tr>
	<td> ct2dot </td>
	<td> ctList to dotbracket </td>
</tr>
<tr>
	<td> write_ctFn </td>
	<td> Save dot-bracket structure to .ct file </td>
</tr>
<tr>
	<td> dot2align </td>
	<td> Convert secondary structure to aligned sequence. </td>
</tr>
<tr>
	<td> dot_from_ctFile </td>
	<td> Read a dotbracket from .ct file </td>
</tr>
<tr>
	<td> trim_stem </td>
	<td> Trim a stem loop </td>
</tr>
<tr>
	<td> find_stem_loop </td>
	<td> Find stem loop from secondary structure </td>
</tr>
<tr>
	<td> find_bulge_interiorLoop </td>
	<td> Find bulges and interior loops from secondary structure </td>
</tr>
<tr>
	<td> calcSHAPEStructureScore </td>
	<td> Calculate strcuture - SHAPE agreement score for stem loop </td>
</tr>
<tr>
	<td> sliding_score_stemloop </td>
	<td> Find stem-loops in RNA with a sliding window </td>
</tr>
<tr>
	<td> multi_alignment </td>
	<td> Multiple sequence alignment with muscle </td>
</tr>
<tr>
	<td> kalign_alignment </td>
	<td> Multiple sequence alignment with kalign </td>
</tr>
<tr>
	<td> global_search </td>
	<td> Global align short sequences to multiple long sequences </td>
</tr>
<tr>
	<td> align_find </td>
	<td> Find the unaligned sequence region from aligned sequence </td>
</tr>
<tr>
	<td> locate_homoseq </td>
	<td> Locate homologous region in multiple sequences </td>
</tr>
<tr>
	<td> dot_to_alignDot </td>
	<td> Dotbracket to aligned dotbracket </td>
</tr>
<tr>
	<td> shape_to_alignSHAPE </td>
	<td> SHAPE list to aligned SHAPE list </td>
</tr>
<tr>
	<td> annotate_covariation </td>
	<td> Annotate raw sequence to colorful sequence by highlight the covariation sites </td>
</tr>
<tr>
	<td> dot_F1 </td>
	<td> Compare predicted structure and true structure and calculate the F1 score </td>
</tr>
<tr>
	<td> parse_structure </td>
	<td> Given a dot-bracket structure, parse structure into all kinds of single-stranded bases and paired bases </td>
</tr>
<tr>
	<td> refine_structure_interior </td>
	<td> Check and make some some canonical base pairs in interior loops paired </td>
</tr>
<tr>
	<td> refine_structure_stackingclosing </td>
	<td> Check and make some some canonical base pairs in stacking end paired </td>
</tr>
<tr>
	<td> refine_structure_hairpinclosing </td>
	<td> Check and make some some canonical base pairs in hairpin paired </td>
</tr>
</table>


<h3> Visual module </h3>
`import Visual`

Prerequisites: java, VARNA (http://varna.lri.fr)

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> Plot_RNAStructure_Shape </td>
	<td> Plot the RNA structure combine with SHAPE scores </td>
</tr>
<tr>
	<td> Plot_RNAStructure_Base </td>
	<td> Plot the RNA structure with different colors for ATCG </td>
</tr>
<tr>
	<td> Plot_RNAStructure_highlight </td>
	<td> Plot the RNA structure and highlight some regions </td>
</tr>
<tr>
	<td> Map_rRNA_Shape </td>
	<td> Output rRNA structure with PostScript format </td>
</tr>
<tr>
	<td> get_rRNA_refseq </td>
	<td> Return reference rRNA sequence </td>
</tr>
</table>

<h3> Rosetta module </h3> 
`from D3 import Rosetta `

Prerequisites: ROSETTA, it can be only run in cluster

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> pred_3D_rosetta </td>
	<td> Predict RNA 3D structure with ROSETTA </td>
</tr>
</table>

<h3> MCSym module </h3> 
`from D3 import MCSym `

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> upload_MCSym_job </td>
	<td> Upload MCSym RNA 3D structure prediction job </td>
</tr>
<tr>
	<td> get_MCSym_status </td>
	<td> Get the status of the job </td>
</tr>
<tr>
	<td> minimize_MCSym_newThread </td>
	<td> Minimize the pdbs </td>
</tr>
<tr>
	<td> score_MCSym_newThread </td>
	<td> Score and ranking pdbs </td>
</tr>
<tr>
	<td> fetch_top_MCSym_pdb_newThread </td>
	<td> Download top scored pdbs </td>
</tr>
</table>

<h3> HDOCK module </h3> 
`from D3 import HDOCK`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> upload_HDOCK_job </td>
	<td> Upload HDOCK RNA-protein docking job </td>
</tr>
<tr>
	<td> get_HDOCK_status </td>
	<td> Get the status of the job </td>
</tr>
<tr>
	<td> guess_HDOCK_time_left </td>
	<td> guess the time to leave </td>
</tr>
<tr>
	<td> fetch_HDOCK_results </td>
	<td> Download all results </td>
</tr>
<tr>
	<td> fetch_HDOCK_top10_results </td>
	<td> Download top 10 results </td>
</tr>
</table>

<h3> Figures module </h3>
`import Figures`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> stackedBarPlot </td>
	<td> Plot a stacked bar figure </td>
</tr>
<tr>
	<td> violinPlot </td>
	<td> Plot a violin figure </td>
</tr>
<tr>
	<td> piePlot </td>
	<td> Plot a pie figure </td>
</tr>
<tr>
	<td> boxPlot </td>
	<td> Plot a box figure </td>
</tr>
<tr>
	<td> cdf </td>
	<td> Plot a CDF curve </td>
</tr>
</table>


### GPU module

`import GPU`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> get_gpu_processes </td>
	<td> Get process handles running on GPU  </td>
</tr>
<tr>
	<td> get_gpu_list </td>
	<td> Get a list of available gpu </td>
</tr>
<tr>
	<td> get_free_gpus </td>
	<td> Get a list of GPU id without process run on it </td>
</tr>
</table>

### Alignment

`import Alignment`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> blast_seq </td>
	<td> Use blastn to search sequence </td>
</tr>
<tr>
	<td> annotate_seq </td>
	<td> Given a sequence and blastdb, search and annotate the sequence </td>
</tr>
</table>

### Covariation

`import Covariation`

<table width="100%">
<tr>
	<th width="20%"> Function name </th>
	<th> Usage </th>
</tr>
<tr>
	<td> dot2sto </td>
	<td> Covert dot to stockholm file </td>
</tr>
<tr>
	<td> cmbuild </td>
	<td> Create .cm file with stockholm alignment </td>
</tr>
<tr>
	<td> cmcalibrate </td>
	<td> Calibrate a .cm file </td>
</tr>
<tr>
	<td> cmsearch </td>
	<td> Call cmsearch programe to search aligned sequence agaist cm model </td>
</tr>
<tr>
	<td> R_scape </td>
	<td> Call R-scape to call covariation base pairs </td>
</tr>
<tr>
	<td> read_RScape_result </td>
	<td> Read the R-scape result </td>
</tr>
<tr>
	<td> get_alignedPos2cleanPos_dict </td>
	<td> Get a distionary {align_pos: raw_pos} </td>
</tr>
<tr>
	<td> call_covariation </td>
	<td> Give sequence and dot. Run covariation pipeline </td>
</tr>
<tr>
	<td> calc_MI </td>
	<td> Calculate the Mutual information for aligned sequences </td>
</tr>
<tr>
	<td> calc_RNAalignfold </td>
	<td> Calculate the RNAalifold covariation score for aligned sequences </td>
</tr>
<tr>
	<td> calc_RNAalignfold_stack </td>
	<td> Calculate the RNAalifold covariation score (consider stack) for aligned sequences </td>
</tr>
<tr>
	<td> collect_columns </td>
	<td> Given multialignment, return alignment columns </td>
</tr>
<tr>
	<td> calc_covBP_from_sto </td>
	<td> Given multialignment, return covariation score for each column pair </td>
</tr>
</table>
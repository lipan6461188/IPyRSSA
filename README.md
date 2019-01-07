## IPyRSSA

IPyRSSA (<b>I</b>ntegrative <b>Py</b>thon library for <b>R</b>NA <b>S</b>econdary <b>S</b>tructure <b>A</b>nalysis) is a set of Python library to analyze RNA secondary structure and SHAPE data.

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
	<td> load_shape </td>
	<td> Read SHAPE .out file </td>
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
</table>

<h3> Structure </h3>

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
	<td> dot2ct </td>
	<td> Dotbracket to list </td>
</tr>
<tr>
	<td> dot2bpmap </td>
	<td> Dotbracket to dictionary </td>
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
</table>










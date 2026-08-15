use clap::{Args, Parser, Subcommand};
// pub(crate) use crate::constants::*;

#[derive(Parser)]
#[clap(author, version, about = "Coverage-adjusted ANI queries and taxonomic profiling for 2bRAD and shotgun metagenomic data.\n\n--- Preparing inputs by extracting (indexing) 2bRAD tags\n## genomes and reads together; outputs <name>.syldb / <name>.sylsp in the -d directory\nmeta2bseek extract -t 5 -g genome1.fa genome2.fa -r sample1.fq -o out -d out -n run1 -e BcgI\n\n## multi-enzyme\nmeta2bseek extract -t 5 -g genome1.fa -r sample1.fq -o out -d out -n run1 -e BcgI,BslFI,CspCI,AloI\n\n--- Nearest neighbour coverage-adjusted ANI\nmeta2bseek query --db-file out/run1.syldb --sample-file out/run1.sylsp > query.tsv\n\n--- Taxonomic profiling with relative abundances and ANI\nmeta2bseek profile --db-file out/run1.syldb --sample-file out/run1.sylsp --minimum-ani 90 --tsv-name abund.tsv --threads 5 > profile.tsv", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// extract 2bRAD tags from genomes and reads. Genome set -> <name>.syldb; each sample -> <name>.sylsp. 
    #[clap(display_order = 1)]
    Extract(ExtractArgs),
    /// sketch sequences using k-mer sampling (similar to sylph). Each sample.fq -> sample.sp. All *.fa -> *.db.
    #[clap(display_order = 2)]
    Sketch(SketchArgs),
    /// Coverage-adjusted ANI querying between databases and samples.
    #[clap(display_order = 4)]
    Query(ContainArgs),
    ///Species-level taxonomic profiling with abundances and ANIs. 
    #[clap(display_order = 3)]
    Profile(ProfileArgs),
    ///Inspect extracted .db and .sp files.
    #[clap(arg_required_else_help = true, display_order = 5)]
    Inspect(InspectArgs),
    ///View sketched .db and .sp files (Meta2bseek sketch format).
    #[clap(arg_required_else_help = true, display_order = 6)]
    View(ViewArgs),
    ///Mark unique (taxa-specific) tags in .db files.
    #[clap(arg_required_else_help = true, display_order = 7)]
    Mark(MarkArgs),
    ///Merge multiple .syldb database files into one .syldb (e.g. per-chunk batch builds).
    #[clap(arg_required_else_help = true, display_order = 8)]
    Merge(MergeArgs),
}


#[derive(Args, Default)]
pub struct ExtractArgs {
    #[clap(short='g', long="genomes", num_args=1.., help_heading = "GENOME INPUT", help = "One or more genome files in fasta format")]
    pub genomes: Option<Vec<String>>,

    #[clap(short='k', long="genome-list", help_heading = "GENOME INPUT", help = "Genome list file: two-column TSV (genome_id<TAB>fasta_path) or one path per line")]
    pub genome_list: Option<String>,

    #[clap(short='r', long="reads", num_args=1.., help_heading = "READ INPUT", help = "One or more fastq files for reads")]
    pub reads: Option<Vec<String>>,

    #[clap(short='s', long="sample-list", help_heading = "READ INPUT", help = "Text file containing paths to fastq files (one per line)")]
    pub sample_list: Option<String>,

    #[clap(short='o', long="output", default_value = ".", help_heading = "OUTPUT", help = "Output directory for extracted tags")]
    pub output_dir: String,

    #[clap(short='e', long="enzyme", default_value = "BcgI", help_heading = "ALGORITHM", help = "Restriction enzyme to use")]
    pub enzyme: String,

    #[clap(short='t', long="threads", default_value_t = 3, help = "Number of threads")]
    pub threads: usize,

    #[clap(short='f', long="format", default_value = "fa", help = "Output format (fa or fq)")]
    pub format: String,

    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(short='1', long="first-pair", help_heading = "PAIRED READ INPUT", help = "First pair of paired-end reads")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second-pair", help_heading = "PAIRED READ INPUT", help = "Second pair of paired-end reads")]
    pub second_pair: Vec<String>,

    #[clap(short='d', long="sample-output-dir", help_heading = "OUTPUT", help = "Output directory for sample files")]
    pub sample_output_dir: String,

    #[clap(short='n', long="out-name", help_heading = "OUTPUT", help = "Output name for generated files")]
    pub out_name: Option<String>,

    #[clap(long="l1", help_heading = "BATCH PAIRED READ INPUT", help = "Text file containing paths to first pair of paired-end reads (one per line)")]
    pub first_pair_list: Option<String>,

    #[clap(long="l2", help_heading = "BATCH PAIRED READ INPUT", help = "Text file containing paths to second pair of paired-end reads (one per line)")]
    pub second_pair_list: Option<String>,

    #[clap(short='l', long="list-sequence", help_heading = "INPUT", help = "File containing list of input sequences")]
    pub list_sequence: Option<String>,

    #[clap(long="max-ram", help_heading = "MEMORY", help = "Maximum RAM usage in GB (default: 75% of total system memory)")]
    pub max_ram: Option<usize>,

    #[clap(long="no-tag-seqs", help_heading = "ALGORITHM", help = "Do not store tag sequences in the .syldb (saves ~4-5x memory/disk for large databases; --mismatch 1 will be unavailable)")]
    pub no_tag_seqs: bool,
}

#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(short='g', long="genomes", num_args=1.., help_heading = "GENOME INPUT", help = "One or more genome files in fasta format")]
    pub genomes: Option<Vec<String>>,

    #[clap(short='k', long="genome-list", help_heading = "GENOME INPUT", help = "Genome list file: two-column TSV (genome_id<TAB>fasta_path) or one path per line")]
    pub genome_list: Option<String>,

    #[clap(short='r', long="reads", num_args=1.., help_heading = "READ INPUT", help = "One or more fastq files for reads")]
    pub reads: Option<Vec<String>>,

    #[clap(short='s', long="sample-list", help_heading = "READ INPUT", help = "Text file containing paths to fastq files (one per line)")]
    pub sample_list: Option<String>,

    #[clap(short='1', long="first", num_args=1.., help_heading = "PAIRED READ INPUT", help = "First pair of paired-end reads")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second", num_args=1.., help_heading = "PAIRED READ INPUT", help = "Second pair of paired-end reads")]
    pub second_pair: Vec<String>,

    #[clap(long="l1", help_heading = "BATCH PAIRED READ INPUT", help = "Text file containing paths to first pair of paired-end reads (one per line)")]
    pub first_pair_list: Option<String>,

    #[clap(long="l2", help_heading = "BATCH PAIRED READ INPUT", help = "Text file containing paths to second pair of paired-end reads (one per line)")]
    pub second_pair_list: Option<String>,

    #[clap(short='o', long="output", default_value = ".", help_heading = "OUTPUT", help = "Output directory for genome sketches (.db files)")]
    pub output_dir: String,

    #[clap(short='d', long="sample-dir", default_value = ".", help_heading = "OUTPUT", help = "Output directory for sample sketches (.sp files)")]
    pub sample_output_dir: String,

    #[clap(short='t', long="threads", default_value_t = 3, help = "Number of threads")]
    pub threads: usize,

    #[clap(short='c', long="c-value", default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,

    #[clap(long="k-size", default_value_t = 31, help_heading = "ALGORITHM", help = "K-mer size")]
    pub k: usize,

    #[clap(long="min-spacing", help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers (default: 30; forced to 0 by --sylph-faithful)")]
    pub min_spacing_kmer: Option<usize>,

    #[clap(long="sylph-faithful", help_heading = "ALGORITHM", help = "Build sylph-equivalent sketches for benchmarking against sylph: keep one copy of k-mers repeated in a genome (dedup instead of deleting them) and force --min-spacing 0. Cannot be combined with an explicit non-zero --min-spacing.")]
    pub sylph_faithful: bool,

    #[clap(long="individual", help_heading = "ALGORITHM", help = "Sketch each contig individually")]
    pub individual: bool,

    #[clap(long="no-dedup", help_heading = "ALGORITHM", help = "Disable deduplication")]
    pub no_dedup: bool,

    #[clap(long="fpr", default_value_t = 0.001, help_heading = "ALGORITHM", help = "False positive rate for deduplication")]
    pub fpr: f64,

    #[clap(long="no-pseudotax", help_heading = "ALGORITHM", help = "Disable pseudotaxonomy tracking")]
    pub no_pseudotax: bool,

    #[clap(long="max-ram", help_heading = "MEMORY", help = "Maximum RAM usage in GB")]
    pub max_ram: Option<usize>,

    #[clap(long="sample-names", num_args=1.., help_heading = "NAMING", help = "Sample names for output files")]
    pub sample_names: Option<Vec<String>>,

    #[clap(long="list-sample-names", help_heading = "NAMING", help = "File containing sample names (one per line)")]
    pub list_sample_names: Option<String>,

    #[clap(long="debug", help_heading = "DEBUG", help = "Enable debug logging")]
    pub debug: bool,

    #[clap(long="trace", help_heading = "DEBUG", help = "Enable trace logging")]
    pub trace: bool,

    // 新增参数，用于指定数据库输出文件名
    #[clap(long="db-out-name", default_value = "database", help_heading = "OUTPUT", help = "Output database file name")]
    pub db_out_name: String,

    // 指定合并文件的名称
    #[clap(long="out-name", help_heading = "OUTPUT", help = "Name for merged output files (default: 'merged_database' for genomes, 'merged_samples' for reads)")]
    pub out_name: Option<String>,

    // 占位符，与extract保持一致的接口
    #[clap(short='l', long="list-sequence", help_heading = "INPUT", help = "File containing list of input sequences")]
    pub list_sequence: Option<String>,

    // 用于兼容性的字段
    pub files: Vec<String>,
}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(num_args=1.., help = "Pre-extracted *.db/*.sp files. Raw single-end fastq/fasta are allowed and will be automatically extracted to .sp/.db")]
    pub files: Vec<String>,

    #[clap(short='l', long="list", help = "Newline delimited file of file inputs",help_heading = "INPUT/OUTPUT")]
    pub file_list: Option<String>,

    #[clap(long, default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum 2bRAD tag multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(short='M',long, default_value_t = 50., help_heading = "ALGORITHM", help = "Exclude genomes with less than this number of extracted 2bRAD tags")]
    pub min_number_kmers: f64,
    #[clap(long, default_value_t = crate::contain::MIN_TAGS_FOR_GENOME, help_heading = "ALGORITHM", help = "Exclude genomes with fewer than N matched tags (default 50; lower values increase sensitivity on large dense databases)")]
    pub min_tags_genome: usize,
    #[clap(short, long="minimum-ani", help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 90 for query and 95 for profile. Smaller than 95 for profile will give inaccurate results." )]
    pub minimum_ani: Option<f64>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='s', long="sample-threads", help = "Number of samples to be processed concurrently. Default: (# of total threads / 3) + 1 for profile, 1 for query")]
    pub sample_threads: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(short='u', long="estimate-unknown", help_heading = "ALGORITHM", help = "Estimate true coverage and scale sequence abundance in `profile` by estimated unknown sequence percentage" )]
    pub estimate_unknown: bool,

    #[clap(long="mismatch", default_value_t = 0, help_heading = "ALGORITHM", help = "Maximum Hamming mismatches allowed when matching a reference 2bRAD tag to sample tags (0 or 1). 1 increases sensitivity to sequencing errors/strain variation but requires tag sequences in the database.")]
    pub mismatch: usize,

    #[clap(long="sketch-mode", help_heading = "ALGORITHM", help = "Interpret .db/.sp inputs as k-mer sketches (Vec<GenomeSketch> / SequencesSketch) produced by `sketch` instead of 2bRAD tag files, and query them through the same containment/ANI engine. Requires --mismatch 0 (sketches store no tag sequences).")]
    pub sketch_mode: bool,

    #[clap(long="read-error-rate", help_heading = "ALGORITHM", help = "Fixed per-base read error rate used for error-aware ANI inversion in --mismatch 1 mode. Overrides automatic estimation from top candidate genomes.")]
    pub read_error_rate: Option<f64>,

    #[clap(long="no-error-correction", help_heading = "ALGORITHM", help = "Disable error-aware ANI inversion in --mismatch 1 mode (restores the old divergence-only inversion).")]
    pub no_error_correction: bool,

    #[clap(long="allow-enzyme-mismatch", help_heading = "ALGORITHM", help = "Proceed even if the sample sketch was extracted with a different enzyme set than the database (default: hard error when the sample lacks a DB enzyme).")]
    pub allow_enzyme_mismatch: bool,

    
    #[clap(short='I',long="read-seq-id", help_heading = "ALGORITHM", help = "Sequence identity (%) of reads. Only used in -u option and overrides automatic detection. ")]
    pub seq_id: Option<f64>,

    //#[clap(short='l', long="read-length", help_heading = "ALGORITHM", help = "Read length (single-end length for pairs). Only necessary for short-read coverages when using --estimate-unknown. Not needed for long-reads" )]
    //pub read_length: Option<usize>,

    #[clap(short='R', long="redundancy-threshold", help_heading = "ALGORITHM", help = "Removes redundant genomes up to a rough ANI percentile when profiling", default_value_t = 99.0, hide=true)]
    pub redundant_ani: f64,

    #[clap(short='r',long="reads", num_args=1.., help = "Single-end raw reads (fastx/gzip)", display_order = 1, help_heading = "extracting")]
    pub reads: Vec<String>,

    #[clap(short='1', long="first-pairs", num_args=1.., help = "First pairs for raw paired-end reads (fastx/gzip)", help_heading = "extracting")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second-pairs", num_args=1.., help = "Second pairs for raw paired-end reads (fastx/gzip)", help_heading = "extracting")]
    pub second_pair: Vec<String>,

    #[clap(short, default_value_t = 200, help_heading = "extracting", help = "Subsampling rate. Does nothing for pre-extracted files")]
    pub c: usize,
    #[clap(short,long="individual-records", help_heading = "extracting", help = "Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-extracted files")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "extracting", help = "Minimum spacing between selected 2bRAD tags on the database genomes. Does nothing for pre-extracted files")]
    pub min_spacing_kmer: usize,

    #[clap(short='o',long="output-file", help = "Output to this file (TSV format). [default: stdout]", help_heading="INPUT/OUTPUT")]
    pub out_file_name: Option<String>,
    #[clap(long="log-reassignments", help = "Output information for how 2bRAD tags for genomes are reassigned during `profile`. Caution: can be verbose and slows down computation.")]
    pub log_reassignments: bool,


    //Hidden options that are embedded in the args but no longer used... 
    #[clap(short, hide=true, long="pseudotax", help_heading = "ALGORITHM", help = "Pseudo taxonomic classification mode. This removes shared 2bRAD tags between species by assigning 2bRAD tags to the highest ANI species. Requires extractes with --enable-pseudotax option" )]
    pub pseudotax: bool,
    #[clap(long="ratio", hide=true)]
    pub ratio: bool,
    #[clap(long="mme", hide=true)]
    pub mme: bool,
    #[clap(long="mle", hide=true)]
    pub mle: bool,
    #[clap(long="nb", hide=true)]
    pub nb: bool,
    #[clap(long="no-ci", help = "Do not output confidence intervals", hide=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hide=true)]
    pub no_adj: bool,
    #[clap(long="mean-coverage", help_heading = "ALGORITHM", help = "Use the robust mean coverage estimator instead of median estimator", hide=true )]
    pub mean_coverage: bool,

}

#[derive(Args)]
pub struct InspectArgs {
    #[clap(num_args=1.., help = "Pre-extracted *.db/*.sp files.")]
    pub files: Vec<String>,
    #[clap(short='o',long="output-file", help = "Output to this file (YAML format). [default: stdout]")]
    pub out_file_name: Option<String>,
    #[clap(long="log-path", help = "Path to store TSV output files")]
    pub log_path: Option<String>,
    #[clap(long="tsv-name", default_value = "tag_matrix.tsv", help = "Name of the TSV file for tag count matrix")]
    pub tsv_name: String,
}

#[derive(Parser, Debug)]
pub struct ProfileArgs {
    #[arg(long)]
    pub sample_file: String,
    
    #[arg(long)]
    pub db_file: String,
    
    #[arg(long, help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 95 for profile. Smaller than 95 for profile will give inaccurate results.")]
    pub minimum_ani: Option<f64>,
    
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    #[arg(long, default_value_t = 0, help_heading = "ALGORITHM", help = "Maximum Hamming mismatches allowed when matching a reference 2bRAD tag to sample tags (0 or 1).")]
    pub mismatch: usize,

    #[arg(long, help_heading = "ALGORITHM", help = "Interpret --db-file as a k-mer sketch database (.db: Vec<GenomeSketch>) and --sample-file as sample sketches (.sp: SequencesSketch or merged Vec<SequencesSketch>) produced by `sketch`, and profile them through the same containment/ANI engine. Requires --mismatch 0 (sketches store no tag sequences).")]
    pub sketch_mode: bool,

    #[arg(long, help_heading = "ALGORITHM", help = "Fixed per-base read error rate used for error-aware ANI inversion in --mismatch 1 mode. Overrides automatic estimation from top candidate genomes.")]
    pub read_error_rate: Option<f64>,

    #[arg(long, help_heading = "ALGORITHM", help = "Disable error-aware ANI inversion in --mismatch 1 mode (restores the old divergence-only inversion).")]
    pub no_error_correction: bool,

    #[arg(long, help_heading = "ALGORITHM", help = "Proceed even if the sample sketch was extracted with a different enzyme set than the database (default: hard error when the sample lacks a DB enzyme).")]
    pub allow_enzyme_mismatch: bool,
    
    #[arg(long)]
    pub out_file_name: Option<String>,

    #[arg(long)]
    pub log_path: Option<String>,

    #[arg(long, default_value = "abundance_matrix.tsv")]
    pub tsv_name: String,
    
    #[arg(long, help = "Taxonomy annotation file (e.g., taxonomy.txt) for species-level aggregation")]
    pub taxonomy_file: Option<String>,
    
    #[arg(long, default_value_t = 10.0, help_heading = "ALGORITHM", help = "Minimum G-score threshold for species filtering. G-score = sqrt(reads_count * tag_count). Default is 10.0")]
    pub gscore_threshold: f64,
    
    #[arg(long, default_value_t = 0.0, help_heading = "ALGORITHM", help = "Minimum effective coverage (eff_cov) for a genome to be retained in profile output. Default 0.0 = no filter.")]
    pub min_eff_coverage: f64,

    #[arg(long, default_value_t = crate::contain::MIN_SHARED_TAGS_FOR_PROFILE, help_heading = "ALGORITHM", help = "Minimum shared tags required to report a genome in profile (default 100; lower values increase sensitivity on large dense databases)")]
    pub min_shared_tags: usize,

    #[arg(long, default_value_t = crate::contain::MIN_TAGS_FOR_GENOME, help_heading = "ALGORITHM", help = "Exclude genomes with fewer than N matched tags (default 50; lower values increase sensitivity on large dense databases)")]
    pub min_tags_genome: usize,

    #[arg(
        long,
        help_heading = "ALGORITHM",
        help = "Enable winner-take-all over-stripping protection (default: OFF). When ON, a genome with >=100 initially shared tags that keeps <=5% of them after reassignment (near-total annihilation) and would fall below the reporting floor is kept with its pre-reassignment ANI/tag counts; its abundance still uses post-reassignment coverage so reads are not double-counted. NOTE: on large dense databases this rescues many strain-cluster losers alongside true positives (precision cost measured on CAMI2); OFF is the safer default."
    )]
    pub reassign_protection: bool,
}

#[derive(Debug)]
pub struct QueryArgs {
    pub query_file: String,
    pub db_file: String,
    pub out_file_name: Option<String>,
    pub minimum_ani: f64,
    pub threads: usize,
}

#[derive(Args)]
pub struct ViewArgs {
    #[clap(num_args=1.., help = "Pre-sketched *.db/*.sp files (Meta2bseek sketch format).")]
    pub files: Vec<String>,
    #[clap(short='o',long="output-file", help = "Output to this file (YAML format). [default: stdout]")]
    pub out_file_name: Option<String>,
    #[clap(long="log-path", help = "Path to store TSV output files")]
    pub log_path: Option<String>,
    #[clap(long="tsv-name", default_value = "kmer_matrix.tsv", help = "Name of the TSV file for k-mer count matrix")]
    pub tsv_name: String,
}

#[derive(Args)]
pub struct MergeArgs {
    #[clap(num_args=1.., required = true, help = "Input .syldb files to merge (must be built with the same enzyme set)")]
    pub inputs: Vec<String>,

    #[clap(short='o', long="output", required = true, help_heading = "OUTPUT", help = "Output merged .syldb file")]
    pub output: String,
}

#[derive(Args)]
pub struct MarkArgs {
    #[clap(help = "Input .db file to mark unique tags")]
    pub input_file: String,
    
    #[clap(short='o', long="output", help = "Output .db file with unique tags marked. If not specified, overwrites input file")]
    pub output_file: Option<String>,
    
    #[clap(long="debug", help = "Enable debug output")]
    pub debug: bool,
}
    
import os
from selina.preprocess.python.scRNA_qc import scrna_qc


def preprocess_parser(subparsers):
    """
    Add argument parsers.
    """

    workflow = subparsers.add_parser(
        "preprocess",
        help="Run preprocess pipeline from scRNA-seq gene-cell count matrix.")
    group_input = workflow.add_argument_group("Input files arguments")
    group_input.add_argument("--format",
                             dest="format",
                             default="",
                             choices=["h5", "mtx", "plain"],
                             help="Format of the count matrix file.")
    group_input.add_argument(
        "--matrix",
        dest="matrix",
        default="",
        help="Location of count matrix file. "
        "If the format is 'h5' or 'plain', users need to specify the name of the count matrix file."
        "If the format is 'mtx', the 'matrix' should be the name of .mtx formatted matrix file, such as 'matrix.mtx'."
    )
    group_input.add_argument(
        "--separator",
        dest="separator",
        default="tab",
        choices=["tab", "space", "comma"],
        help="The separating character (only for the format of 'plain')."
        "Values on each line of the plain matrix file will be separated by the character. DEFAULT: tab."
    )
    group_input.add_argument(
        "--feature",
        dest="feature",
        default="features.tsv",
        help="Location of feature file (required for the format of 'mtx'). "
        "Features correspond to row indices of count matrix. DEFAULT: features.tsv."
    )
    group_input.add_argument(
        "--gene-column",
        dest="gene_column",
        default=2,
        type=int,
        help=
        "If the format is 'mtx', please specify which column of the feature file to use for gene names. DEFAULT: 2."
    )
    group_input.add_argument(
        "--gene-idtype",
        dest="gene_idtype",
        default="symbol",
        choices=["symbol", "ensembl"],
        help=
        "Type of gene name, 'symbol' for gene symbol and 'ensembl' for ensembl id. DEFAULT: symbol."
    )
    group_input.add_argument(
        "--barcode",
        dest="barcode",
        default="barcodes.tsv",
        help="Location of barcode file (required for the format of 'mtx'). "
        "Cell barcodes correspond to column indices of count matrix. DEFAULT: barcodes.tsv. "
    )
    # group_input.add_argument(
    #     "--meta-file",
    #     dest="meta_file",
    #     default="",
    #     help="Location of metadata file. "
    #     "The metadata file should be a table with cells as rows and meta information as columns. "
    #     "The first line of the metadata file should contain the names of the variables."
    # )
    # group_input.add_argument(
    #     "--meta-sep",
    #     dest="meta_sep",
    #     default="tab",
    #     choices=["tab", "space", "comma"],
    #     help="The separating character of the metadata file. "
    #     "Values on each line of the metadata file will be separated by the character. DEFAULT: tab."
    # )
    # group_input.add_argument(
    #     "--meta-cell-column",
    #     dest="meta_cell",
    #     default=1,
    #     type=int,
    #     help=
    #     "Please specify which column of the metadata file to use for cell ID. DEFAULT: 1."
    # )
    group_input.add_argument(
        "--assembly",
        dest="assembly",
        default="GRCh38",
        choices=["GRCh38", "GRCh37"],
        type=str,
        help="Assembly (GRCh38/hg38 and GRCh37/hg19). DEFAULT: GRCh38.")

    # Quality control cutoff
    group_cutoff = workflow.add_argument_group("Quality control arguments")
    group_cutoff.add_argument(
        "--count-cutoff",
        dest="count_cutoff",
        default=1000,
        type=int,
        help="Cutoff for the number of count in each cell. DEFAULT: 1000.")
    group_cutoff.add_argument(
        "--gene-cutoff",
        dest="gene_cutoff",
        default=500,
        type=int,
        help=
        "Cutoff for the number of genes included in each cell. DEFAULT: 500.")
    group_cutoff.add_argument(
        "--cell-cutoff",
        dest="cell_cutoff",
        default=10,
        type=int,
        help="Cutoff for the number of cells covered by each gene. DEFAULT: 10."
    )

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument(
        "--directory",
        dest="directory",
        default="preoprocess",
        help=
        "Path to the directory where the result file shall be stored. DEFAULT: preprocess."
    )
    group_output.add_argument("--outprefix",
                              dest="outprefix",
                              default="query",
                              help="Prefix of output files. DEFAULT: query.")
    group_output.add_argument(
        "--mode",
        dest="mode",
        choices=["single", "cluster", 'both'],
        required=True,
        help=
        "Output expression file for prediction. single: single-cell level. cluster: cluster level. both: output both the single-cell level and cluster level expression profiles"
    )


# Generate Rscript
# def GenerateRscript(count_file, gene_idtype, gene_cutoff, cell_cutoff,
#                     meta_file, meta_sep, meta_cell, assembly, outprefix,
#                     directory):


# Generate Rscript
def GenerateRscript(count_file, gene_idtype, gene_cutoff, cell_cutoff,
                    assembly, outprefix, directory, mode):

    rfile = os.path.join(directory, "%s.R" % (outprefix))
    outf = open(rfile, "w")
    refGenePath = os.path.split(os.path.split(
        os.path.abspath(__file__))[0])[0] + '/data'
    rsrcPath = os.path.split(os.path.split(
        os.path.abspath(__file__))[0])[0] + '/r'

    # load package
    script = '''# load package
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
source("%s/RNARunSeurat.R")
source('%s/RNAAssemblyConvert.R')
source('%s/RNAEnsemblToSymbol.R')
source('%s/FindAllMarkersMAESTRO.R')
\n
''' % (rsrcPath, rsrcPath, rsrcPath, rsrcPath)
    outf.write(script)

    # read data
    script = '''# read data
expr = Read10X_h5("%s")
''' % (count_file)
    outf.write(script)

    # assembly conversion and gene id conversion
    if assembly == "GRCh37":
        if gene_idtype == "symbol":
            script = '''
# assembly conversion
expr = RNAAssemblyConvert(expr, from = "GRCh37", to = "GRCh38", dataPath = "%s")
''' % (refGenePath)
        elif gene_idtype == "ensembl":
            script = '''
# gene id conversion
expr = RNAEnsemblToSymbol(expr, assembly = "GRCh37", dataPath = "%s")
expr = RNAAssemblyConvert(expr, from = "GRCh37", to = "GRCh38", dataPath = "%s")
''' % (refGenePath, refGenePath)
        outf.write(script)

    elif assembly == "GRCh38":
        if gene_idtype == "ensembl":
            script = '''
# gene id conversion
expr = RNAEnsemblToSymbol(expr, assembly = "GRCh38", dataPath = "%s")
''' % (refGenePath)
            outf.write(script)

    assembly = "GRCh38"

    # analysis
    script = '''
# choose optimal pc and resolution based on cell number
cells <- ncol(expr)
if (cells <= 5000) {
npc <- 50; cluster.res <- 0.6
} else if(cells <= 100000) {
npc <- 50; cluster.res <- 1
} else {
npc <- 100; cluster.res <- 1
}\n

# clustering
RNA.res = RNARunSeurat(inputMat = expr, 
                       project = "%s", 
                       min.c = %d,
                       min.g = %d,
                       npcs = npc,
                       variable.genes = 2000, 
                       assembly = "%s",
                       cluster.res = cluster.res,
                       outdir = "%s",
                       mode = "%s")\n
''' % (outprefix, cell_cutoff, gene_cutoff, assembly, directory, mode)
    outf.write(script)

    # # read metadata
    # if meta_file:
    #     if meta_sep == "tab":
    #         sep = "\\t"
    #     elif meta_sep == "space":
    #         sep = " "
    #     elif meta_sep == "comma":
    #         sep = ","
    #     script = '''
    # add metadata
    # meta = read.delim("%s", header = T, row.names = %d, sep = "%s")
    # RNA.res$RNA@meta.data = cbind(RNA.res$RNA@meta.data, meta[colnames(RNA.res$RNA),, drop = FALSE])
    # for (i in colnames(meta)) {
    #   p = DimPlot(object = RNA.res$RNA, group = i, label = FALSE, pt.size = 0.1)
    #   if (length(unique(RNA.res$RNA@meta.data[,i])) > 20) {
    #     plot_width = 10
    #   } else {
    #     plot_width = 7
    #   }
    #   ggsave(file.path("%s", paste0(RNA.res$RNA@project.name, "_", i, ".png")), p,  width=plot_width, height=5)
    # }
    # ''' % (os.path.abspath(meta_file), meta_cell, sep, directory)
    #         outf.write(script)

    script = '''
# save object
saveRDS(RNA.res, "%s_res.rds")
''' % (os.path.join(directory, outprefix))
    outf.write(script)

    outf.close()
    return os.path.abspath(rfile)


# def query_preprocess(directory, outprefix, fileformat, matrix, separator,
#  feature, gene_column, gene_idtype, barcode, meta_file,
#  meta_sep, meta_cell, count_cutoff, gene_cutoff,
#  cell_cutoff, assembly):


def query_preprocess(directory, outprefix, fileformat, matrix, separator,
                     feature, gene_column, gene_idtype, barcode, count_cutoff,
                     gene_cutoff, cell_cutoff, assembly, mode):

    scrna_qc(directory + "/Data", outprefix, fileformat, matrix, separator,
             feature, gene_column, barcode, count_cutoff, gene_cutoff,
             cell_cutoff, assembly)

    count_file = os.path.abspath(
        os.path.join(directory + "/Data",
                     outprefix + "_filtered_gene_count.h5"))

    # rscript = GenerateRscript(count_file, gene_idtype, gene_cutoff,
    #                           cell_cutoff, meta_file, meta_sep, meta_cell,
    #                           assembly, outprefix, directory)
    rscript = GenerateRscript(count_file, gene_idtype, gene_cutoff,
                              cell_cutoff, assembly, outprefix, directory,
                              mode)

    cmd = "Rscript %s" % (rscript)
    os.system(cmd)
    cmd = "rm -r %s; rm -r %s" % (rscript, directory + "/Data")
    os.system(cmd)
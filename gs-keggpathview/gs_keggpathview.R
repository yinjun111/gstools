library("argparser",quietly =T)

version="0.7"

# v0.1: pathview function
# v0.2: TST genelist, Log2FC limit
# v0.3: Argument requirements, pathway lists, auto comparison name, outfile location
# v0.4: Graphviz view (PDF) enabled
# v0.5: Geneanno not required
# v0.6: --defile is now --in, GeneDE input, --fccol, updated name gen, out folder creation
# v0.6: Gstools --gs, --top, KEGG reference file, --showname, --pathwayname search
# v0.7: New KEGG pathway list, updated --pathwayname search and --showname output

description=paste0("gs_keggpathview\nversion ",version,"\n","Usage:\nDescription: Generate KEGG pathway images for DE data.\n")


# --------------
# Arguments
# --------------
parser <- arg_parser(description=description)
parser <- add_argument(parser, arg="--in",short="-i", type="character", help = "Input DE file from omictools summary folder, ex. GeneDEMerged.txt, GeneDEreformated.txt")
parser <- add_argument(parser, arg="--fccol", short="-f", type="numeric", default=2, help = "Log2FC column number to use, counting rownames as a column (Default = 2)")
parser <- add_argument(parser, arg="--name",short="-n", type="character", help = "Custom name for the output file")
parser <- add_argument(parser, arg="--out",short="-o", type="character", help = "Output folder, png format")
parser <- add_argument(parser, arg="--pathwayid",short="-p", type="character",help = "KEGG Pathway ID(s) to be visualized, ex. 04110")
parser <- add_argument(parser, arg="--gs", type="character", help = "Input file from gstools KEGG pathway analysis, ex. c2.cp.kegg/Comparison_gs-fisher.txt")
parser <- add_argument(parser, arg="--top", type="numeric", help = "How many top pathways to plot. Use with --gs")
parser <- add_argument(parser, arg="--pathwayname", type="character", help = "Choose pathway(s) by name, separated by comma, with underscores \'_\' instead of spaces")
parser <- add_argument(parser, arg="--showname", flag=TRUE, help = "Whether or not to include the pathway name in the outfile name (logical).")
#parser <- add_argument(parser, arg="--genes",short="-g", type="character",help = "Genes to be shown, separated by comma ',' (not yet implemented)")
#parser <- add_argument(parser, arg="--geneanno",short="-a", type="character", help = "Ensembl ID and Gene ID file needed for gene list (--genes), ex. geneanno.txt")
parser <- add_argument(parser, arg="--tx",short="-t", type="character",help = "Currently support Human.B38.Ensembl88, Mouse.B38.Ensembl88, Rat.Rn6.Ensembl88")
parser <- add_argument(parser, arg="--fccutoff", type="character", default ="1",help = "Log2FC cutoff to use for mapping color scale. Any abs(Log2FC) values at or above the cutoff will be colored at the max/min on the color scale")
#parser <- add_argument(parser, arg="--qcutoff", type="character", default ="0.05",help = "BH-adjusted p-value significance cutoff. Removes nonsignificant genes")
parser <- add_argument(parser, arg="--outtype", type="character", default = "pdf", help = "Export 'png' or 'pdf' pathway map")

args = parse_args(parser)

# -----------------
# Kegg Reference
# -----------------
# Pcluster Slurm
kegg.path <- "/data/jyin/Databases/gstools-db/KEGG/kegg.pathwaylist.2022.10.28.txt"

# --------------
# Requirements
# --------------
# Infile Requirement (DEA)
  # Should be adjusted when adding genes functionality
if (is.na(args$"in")) {
  stop("Must specify a file for --in (GeneDEMerged)")
}

# Genes not yet implemented
# if (!is.na(args$"genes")) {
#   stop("Gene list not yet implemented for --genes")
# }

# Genelist Warning
# if (!is.na(args$"in")) {
#   if (!is.na(args$"genes")) {
#     stop("Cannot use --in and --genes together")
#   }
# }

# Fccol Check
if (args$"fccol" < 2) {
  stop("Must set --fccol to 2 or greater.")
}
fccol <- as.numeric(args$"fccol" - 1) # Adjusting fccol since column 1 becomes rownames

colcheck <- fccol + 4
if (colcheck %% 5 != 0) {
  print(paste0("Warning: Fold Change column (--fccol) set to ", args$"fccol",". This may not be a Log2FC column. Consider setting to 2, 7, 12, etc."))
}

# Outtype Options
type <- FALSE
if (!is.na(args$"outtype")) {
  if (args$"outtype" == "pdf") {
    type <- FALSE
    print("Graphviz PDF selected as output type.")
  } else if (args$"outtype" == "png") {
    type <- TRUE
    print("KEGG Native png selected as output type.")
  } else {
    stop("Must set --outtype to either pdf or png")
  }
}

# Pathway List
if (!is.na(args$"gs")) {
  if (!is.na(args$"pathwayname")) {
    print("Warning: Using pathways from --gs gs-fisher input. Ignoring --pathwayname pathway list")
  }
  if (!is.na(args$"pathwayid")) {
    print("Warning: Using pathways from --gs gs-fisher input. Ignoring --pathwayid pathway list")
  }
} else if (is.na(args$"gs")) {
  if (!is.na(args$"pathwayname")) {
    if (!is.na(args$"pathwayid")) {
      stop(paste0("Cannot use --pathwayname with --pathwayid, see list of pathway IDs and pathway names at ",kegg.ref.path))
    }
  }
}

# Out Folder Requirement and Creation
if (is.na(args$"out")) {
  print(paste0("Warning: Output folder not specified using --out, generating outputs in local script folder ", getwd()))
} else if (!is.na(args$"out")) {
  # Check if folder exists
  if (file.exists(args$"out")) {
    print(paste0("Generating outputs in ", args$"out"))
  } else if (!file.exists(args$"out")) {
    # Create folder if does not exist
    print(paste0("Creating folder ", args$"out"))
    dir.create(args$"out")
  }
}

# Species Code Requirement
# Kegg Reference Paths
if (is.na(args$"tx")) {
  stop("Need to specify a species transcriptome with --tx")
} else if (args$tx == "Human.B38.Ensembl88") {
  species <- "hsa"
} else if (args$tx == "Mouse.B38.Ensembl88") {
  species <- "mmu"
} else if (args$tx == "Rat.Rn6.Ensembl88") {
  species <- "rno"
} else {
  stop("Species transcriptome not supported, see --help for --tx")
}

# --------------
# Packages
# --------------

library(pathview)
library(dplyr)

# --------------
# Defaults
# --------------
# Colors: Low, Mid, High
cols <- c("#0571b0", "gray", "#ca0020")
# Gene ID Type
id.type <- "ENSEMBL"

# --------------
# Log2FC Cutoff
# --------------
fccutoff <- as.numeric(args$"fccutoff")
if (is.na(fccutoff)) {
  stop("Must specify numeric fccutoff")
}
print(paste0("Using Log2FC cutoff ", fccutoff, " for map color scale min/max"))

# --------------
# Input Data
# --------------
# Data
data <- read.table(args$"in", header=T, row.names=1, sep="\t", quote="", comment.char="")

# Assigning --name or Automatically grabbing comparison name if using GeneDE
  # For getting name from file name instead of column name
    # depath <- unlist(strsplit(args$"in", split = "/"))
    # dename <- depath[length(depath)]
    # desplit <- unlist(strsplit(dename, split = "_"))
    # dename <- paste0(desplit[-length(desplit)], collapse = "_")
if (is.na(args$"name")) {
  if (grepl("GeneDE",args$"in", ignore.case = TRUE)) {
    # Changing Log2FC Group1 vs Group2 to Group1_vs_Group2
    dename <- colnames(data)[fccol]
    dename <- gsub(".vs.","_vs_", dename, ignore.case = TRUE)
    dename <- gsub("Log2FC.", "", dename, ignore.case = TRUE)
    
  } else {
    print("No --name specified, outfile will only contain species and pathway details.")
    dename <- ""   
  }
} else if (!is.na(args$"name")) {
  dename <- args$"name"
}

# Gene IDs not yet implemented
# if (!is.na(args$"genes")) {
#   geneanno <- read.table(args$"geneanno", header=T, row.names=1, sep="\t", quote="", comment.char="")
#   # Function: Combining Data & GeneID
#   make.geneid <- function(data = NULL, geneanno = NULL) {
#     geneanno.red <- geneanno[,c(1,2)] # Reduced gene annotation
#     data.anno <- merge(geneanno.red, data, by = 0) # Annotating data
#     rownames(data.anno) <- data.anno$Row.names # Ensembl as Rownames
#     data.anno$Row.names <- NULL # Formatting
#     data.anno$description <- NULL # Formatting
#     data.anno = data.anno %>% dplyr::select(-gene_name, gene_name)
#     data.anno$gene_name <- toupper(data.anno$gene_name)
#     return(data.anno)
#   }
#   # Combining Data & GeneID
#   data <- make.geneid(data = data, geneanno = geneanno)
# }

# -----------
# Gene List
# -----------
# Not yet implemented

# Subsetting for Gene List
# if (!is.na(args$"genes")) {
#   # Formatting
#   genes <- gsub(" ", "", args$"genes")
#   genes <- strsplit(genes, split = ",")
#   genes <- genes[[1]]
#   genes <- toupper(genes)
#   
#   # Subsetting
#   data <- data[data$gene_name %in% genes,]
#   # Setting all selected genes FC to the cutoff (Red color)
#   data[,1] <- fccutoff
# }

# Green color if present when using gene list
# if (!is.na(args$"genes")) {
#   cols[1] <- "gray"
#   cols[3] <- "#56953d" # green
# }

# Fold Change Vector
data.fc <- as.numeric(as.vector(unlist(data[,fccol])))
names(data.fc) <- rownames(data)

# Significance Vector
#data.sig <- as.numeric(as.vector(unlist(data[,5])))
#names(data.sig) <- rownames(data)

# -------------------
# Gstools Fisher
# Pathway Assignment
# -------------------
# KEGG gs-Annotation to Pathway ID
kegg.ref <- read.table(kegg.path, header=T, row.names=NULL, sep="\t", quote="", comment.char="", colClasses=c("character"))

if (!is.na(args$"gs")) {
  print("Using gs-fisher pathways (--gs).")
  
  # --gs
  gsdata <- read.table(args$"gs", header=T, row.names=NULL, sep="\t", quote="", comment.char="")
  
  # Remove Total Genes row, keep Pathway and Significance columns
  gsdata <- gsdata[-1,]
  gsdata <- gsdata[,c(1,9)]
  
  # Subsetting for Pathways
  if (is.na(args$"top")) {
    gsdata <- gsdata[gsdata[,2] == TRUE,]
    print(paste0("Plotting ",nrow(gsdata)," significant pathway(s) from gs-fisher file."))
  } else if (!is.na(args$"top")) {
    # Top Pathways
    gsdata <- gsdata[c(1:args$"top"),]
    print(paste0("Plotting top ",args$"top"," pathways from gs-fisher file."))
  }
  
  # --gs
  pathways <- kegg.ref[kegg.ref[,1] %in% gsdata[,1],2]
  names(pathways) <- kegg.ref[kegg.ref[,1] %in% gsdata[,1],1]
  pathways <- na.omit(pathways)
  
} else if (is.na(args$"gs")) {
  if (!is.na(args$"pathwayid")) {
    # --pathwayid
    print("Using pathway IDs (--pathwayid).")
    pathways <- unique(unlist(strsplit(args$"pathwayid", split = ",")))
    names(pathways) <- kegg.ref[kegg.ref[,2] %in% pathways,1]
    
  } else if (!is.na(args$"pathwayname")) {
    # --pathwayname
    print("Searching pathway name(s) (--pathwayname).")
    # Old code for MSigDB Pathway Names. New names contain commas, so multiple pathway names cannot be searched at the same time.
    # pathwaynames <- gsub("-","_",args$"pathwayname")
    # pathwaynames <- gsub("'","",pathwaynames)
    # pathwaynames <- gsub("KEGG_","",pathwaynames)
    # pathwaynames <- gsub("/","_",pathwaynames)
    # pathwaynames <- unique(unlist(strsplit(args$"pathwayname", split = ",")))
    # Exact Matches
    # greppaths <- unlist(lapply(args$"pathwayname",grep,kegg.ref[,1],ignore.case = TRUE, value = TRUE))
    # Pathway IDs
    pathways <- kegg.ref[kegg.ref[,1] %in% args$"pathwayname",2]
    # Pathway ID Names
    names(pathways) <- kegg.ref[kegg.ref[,1] %in% args$"pathwayname",1]
    pathways <- na.omit(pathways)
    
    if (length(pathways)==0) {
      stop(paste0("Pathway from --pathwayname does not have a corresponding pathway ID, see list of pathway names at ",kegg.path))
    }
  }
}

# Set Out Folder
if (!is.na(args$"out")) {
  setwd(args$"out")
}
print(paste0("Generating outputs in ", getwd()))

# --------------
# Function
# --------------
for (pathway in pathways) {
  # Adding Path Name to Outfile Name
  if (args$"showname") {
    pathname <- names(pathways[pathways==pathway])
    # Formatting File Name
    pathname <- gsub(" - ","__",pathname)
    pathname <- gsub(" / ","_",pathname)
    pathname <- gsub("-","_",pathname)
    pathname <- gsub(" ","_",pathname)
    pathname <- gsub(",","",pathname)
    pathname <- gsub("\\(","",pathname)
    pathname <- gsub("\\)","",pathname)
    pathname <- substr(pathname, 1,20) # First 20 characters only
    name <- paste0(pathname,".",dename) 
  } else {
    name <- dename
  }
  # Plotting Pathway
  pv.out <- pathview(gene.data = data.fc,
                     gene.idtype = id.type,
                     pathway.id = pathway,
                     species = species,
                     out.suffix = name,
                     kegg.native = type,
                     same.layer = T,
                     low = list(gene = cols[1], cpd = "green"),
                     mid = list(gene = cols[2], cpd = "gray"),
                     high = list(gene = cols[3], cpd = "yellow"),
                     limit = list(gene=fccutoff, cpd=1))
}



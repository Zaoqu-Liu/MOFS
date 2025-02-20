#' @title Perform ssGSEA-based Subtyping Using Established Molecular Markers
#' @description This function performs single-sample Gene Set Enrichment Analysis (ssGSEA) to classify samples into subtypes based on predefined molecular markers, including PN (proneural), MES (mesenchymal), and OXPHOS (oxidative phosphorylation). The function calculates enrichment scores for each subtype and assigns the most likely classification to each sample, facilitating molecular subtype analysis for glioma studies.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param data.test A matrix or data frame representing the input expression data, where rows are genes and columns are samples.
#' @param dir.file Character. Directory for saving the output files (default: '.').
#' @param gct.filename Character. The filename for the generated GCT file (default: 'data.gct').
#' @param number.perms Integer. Number of permutations for ssGSEA analysis (default: 10).
#' @param tolerate.mixed Logical. Whether to allow "Mixed" predictions when multiple gene sets have the same minimum p-value (default: FALSE).
#' @return A data frame with the following columns:
#'   - `ID`: Sample identifiers.
#'   - `Predict`: Predicted subtype for each sample.
#'   - Columns with `_pval`: P-values for each subtype.
#' @details The function uses predefined marker gene sets for GBM subtypes:
#'   - `PN`: Proneural subtype markers.
#'   - `OXPHOS`: Oxidative phosphorylation subtype markers.
#'   - `MES`: Mesenchymal subtype markers.
#'
#' The function integrates `Classification.ssGSEA` to:
#'   1. Perform ssGSEA for the input expression data.
#'   2. Predict subtypes based on the most significant enrichment (smallest p-value).
#'   3. Optionally assign "Mixed" label if multiple subtypes have the same minimum p-value.
#' @references Wu M, Wang T, Ji N, Lu T, Yuan R, Wu L, et al. Multi-omics and pharmacological characterization of patient-derived glioma cell lines. Nat Commun. 2024;15:6740. doi:10.1038/s41467-024-51214-y.
#' @examples
#' # Simulated expression data
#' data.test <- matrix(rnorm(10000), nrow = 100, ncol = 100)
#' rownames(data.test) <- paste0("Gene", 1:100)
#' colnames(data.test) <- paste0("Sample", 1:100)
#'
#' # Run WuGBM subtyping
#' result <- WuGBM(
#'   data.test = data.test,
#'   dir.file = "./results",
#'   gct.filename = "test_data.gct",
#'   number.perms = 10,
#'   tolerate.mixed = TRUE
#' )
#' print(result)
#' @export
WuGBM <- function(
    data.test,
    dir.file = ".",
    gct.filename = "data.gct",
    number.perms = 10,
    tolerate.mixed = FALSE) {
  # Define marker gene sets for GBM subtypes
  geneList <- list(
    PN = c("ZNF266", "NCAM1", "KCND3", "ZNF711", "MBTD1", "SPEN", "DCP2", "PPP1R12B", "ZNF426", "ZNF532", "ZNF562", "NLGN3", "MARK1", "ZFP30", "APC2", "KCNAB3", "ZNF573", "ZNF8", "ZNF124", "POU2F1", "APBB3", "GNAO1", "SGSM3", "ZNF253", "CHD7", "ZNF614", "ZNF175", "ZNF557", "PTPRZ1", "SS18L1", "USP24", "WSCD1", "ZNF337", "ZNF571", "HMBOX1", "PNN", "TMPRSS5", "PURG", "MAP3K1", "SHC2", "LRRTM2", "YPEL1", "SPATA6", "ZNF345", "LRP4", "MAP2", "PHC1", "SP4", "ZNF211", "ZNF221", "CSPG5", "SOX5", "KLF15", "QKI", "ACVR2B", "DHPS", "EFNA2", "FZD3", "ADAM22", "C1orf50", "KCNIP2", "ZFP2", "ZNF493", "BSN", "GAB1", "SCYL3", "TRIM9", "ABAT", "ATF7IP", "DDX17", "ELAVL3", "ELOVL2", "KCNJ10", "STX6", "AGPAT5", "PMP2", "TTYH1", "ZNF230", "ZNF551", "CEP68", "ONECUT2", "PHF21A", "PROX1", "ZNF780B", "MAPT", "RFX3", "ZNF419", "ZNF549", "ARVCF", "CPSF6", "PGAP1", "POGZ", "TRIM24", "ZDHHC11", "ZNF225", "BCAN", "EPB41L5", "LSM14A", "NPAS3", "ZNF227"),
    OXPHOS = c("MRPL2", "TMEM11", "SMARCD2", "BYSL", "CCDC51", "TBL3", "MAD2L1BP", "SCAND1", "MRPL28", "CCT7", "CUTA", "MRPS18B", "SF3B5", "SIVA1", "TUFM", "NUP85", "PPM1G", "C17orf53", "COX5A", "EIF4A3", "FAU", "TXN2", "UQCRC1", "TPI1", "ACD", "CCT3", "LSM2", "MRPS34", "RPUSD2", "SAC3D1", "SNRPC", "AURKB", "EIF3B", "NUDC", "PHB", "CSK", "ATAD3A", "IMP3", "MRPL11", "SMPD2", "AAMP", "CKS1B", "GPS1", "GSTP1", "PSMB3", "RPL35", "SPAG7", "STOML2", "THAP4", "UXT", "AKR7A2", "CXorf56", "FN3KRP", "GFER", "MRPS12", "NFKBIL1", "PPP4C", "SNRPB", "EDC3", "EIF1", "HIGD2A", "MRPS2", "NUBP2", "RPL19", "SNRPD2", "ATIC", "C20orf27", "CDT1", "DPM3", "EIF2B4", "ELAC2", "FBL", "FKBPL", "MRPL22", "MRPS15", "NDUFA8", "NDUFB11", "POP7", "RALY", "RPS21", "RPS5", "THAP7", "TOMM22", "C1orf35", "EIF5A", "GALK1", "HSD17B10", "JMJD6", "NANS", "PES1", "PQBP1", "PRPF4", "PRPF6", "QARS", "RRP9", "THAP11", "TRAPPC2L", "ABT1", "ERAL1", "LSM4"),
    MES = c("SPAG1", "CAST", "MET", "TBK1", "DCBLD2", "EDEM1", "TWSG1", "GEM", "TACC1", "CEP290", "NEDD4", "RAP2C", "WWP1", "CD44", "KCTD9", "CLIP1", "MYH15", "PSMD5", "ERGIC2", "IMPAD1", "RAB27A", "SLC5A3", "GTF2F2", "SAMD9", "PTPN12", "TMED7", "TMF1", "DHX29", "HRH1", "ITGAV", "SYTL2", "TRPC1", "AKAP9", "PARP4", "PDCD6IP", "EIF2AK2", "FEM1C", "GLIPR1", "GOLGA4", "RABL3", "RASSF8", "ANKRD28", "ARNTL2", "BCAP29", "ELOVL5", "GCC2", "ETAA1", "MYBL1", "DIRAS3", "UFM1", "PHEX", "TGOLN2", "TWF1", "CROT", "MTDH", "NR3C1", "SH2D4A", "LAMP2", "PIGK", "SCG2", "CUL4B", "IQGAP1", "PSMD12", "RASA2", "RECQL", "ZMYM5", "CAPZA1", "CAV1", "TBC1D5", "ZMPSTE24", "PICALM", "SEMA3A", "ARL4C", "SYPL1", "UEVLD", "XRCC4", "ENTPD4", "MPHOSPH10", "PPP2R3A", "BRCA2", "HMMR", "MCFD2", "LRRFIP2", "PRKCI", "SH2B3", "KTN1", "ARMC1", "CD2AP", "DNAJC13", "ITGB1", "MSH3", "SLC9A7", "TEAD1", "ZNF654", "ADAM9", "ATF1", "ELOVL6", "KRR1", "MBTPS2", "PSME4")
  )

  # Run ssGSEA-based subtyping
  res <- Classifier.ssGSEA(
    data.test = data.test,
    dir.file = dir.file,
    gct.filename = gct.filename,
    number.perms = number.perms,
    tolerate.mixed = tolerate.mixed,
    marker.list = geneList
  )

  return(res)
}

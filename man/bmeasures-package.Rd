\name{bmeasures-package}
\alias{bmeasures-package}

\docType{package}
\title{
Binary similarity/dissimilarity measures
}
\description{
This package generates the quantities of the OTUs table, calculates the binary similarity/dissimilarity measures between two vectors, and finding the most suitable binary similarity/dissimilarity equations using ROC analysis.
}
\details{
\tabular{ll}{
Package: \tab bmeasures\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-10-13\cr
License: \tab GPL-2\cr
}
\bold{Functions:}\cr
\tabular{ll}{
\code{bmeasures(x, y, method)} \tab Calculates the binary similarity/dissimilarity coefficient between two vectors. \cr
\code{bmeasures_otu(x, y)} \tab Generates the quantities of the Operational Taxonomic Units table. \cr
\code{bmeasures_find(inFile, setSeed=0, numSample=20)} \tab Finding a suitable binary similarity and dissimilarity measures. 
}
}
\author{
Sony H. Wijaya, Farit M. Afendi, Irmanida Batubara, Latifah K. Darusman, Md. Altaf-Ul-Amin, Shigehiko Kanaya\cr
Maintainer: Sony H. Wijaya <\email{sonyhartono@gmail.com}>

}
\references{
[1]. Avcibaş I, Kharrazi M, Memon N, Sankur B: Image steganalysis with binary similarity measures. EURASIP J Appl Signal Processing 2005, 17:2749–2757.\cr
[2]. Baroni-urbani C, Buser MW: Similarity of binary data. Syst Biol 1976, 25:251–259.\cr
[3]. Batagelj V, Bren M: Comparing resemblance measures. J Classif 1995, 12:73–90.\cr
[4]. Boyce RL, Ellison PC: Choosing the best similarity index when performing fuzzy set ordination on binary data. J Veg Sci 2001, 12:711–720.\cr
[5]. Cha S-H, Tappert CC, Yoon S: Enhancing Binary Feature Vector Similarity Measures. 2005.\cr
[6]. Cha S, Choi S, Tappert C: Anomaly between Jaccard and Tanimoto coefficients. In Proceedings of Student-Faculty Research Day, CSIS, Pace University; 2009:1–8.\cr
[7]. Chang J, Chen R, Tsai S: Distance-preserving mappings from binary vectors to permutations. IEEE Trans Inf Theory 2003, 49:1054–1059.\cr
[8]. Cheetham AH, Hazel JE, Journal S, Sep N: Binary (presence-absence) similarity coefficients. J Paleontol 1969, 43:1130–1136.\cr
[9]. Choi S-S, Cha S-H, Tappert CC: A survey of binary similarity and distance measures. J Syst Cybern Informatics 2010, 8:43–48.\cr
[10]. Consonni V, Todeschini R: New similarity coefficients for binary data. Match-Communications Math Comput Chem 2012, 68:581–592.\cr
[11]. da Silva Meyer A, Garcia AAF, Pereira de Souza A, Lopes de Souza C: Comparison of similarity coefficients used for cluster analysis with dominant markers in maize (Zea mays L). Genet Mol Biol 2004, 27:83–91.\cr
[12]. Dalirsefat SB, da Silva Meyer A, Mirhoseini SZ: Comparison of similarity coefficients used for cluster analysis with amplified fragment length polymorphism markers in the silkworm, Bombyx mori. J Insect Sci 2009, 9:1–8.\cr
[13]. Dice LR.: Measures of the amount of ecologic association between species. Ecology 1945, 26:297–302.\cr
[14]. Faith DP: Asymmetric binary similarity measures. Oecologia 1983, 57:287–290.\cr
[15]. Gower JC, Legendre P: Metric and Euclidean properties of dissimilarity coefficients. J Classif 1986, 3:5–48.\cr
[16]. Holliday JD, Hu C-Y, Willett P: Grouping of coefficients for the calculation of inter-molecular similarity and dissimilarity using 2D fragment bit-strings. Comb Chem High Throughput Screen 2002, 5:155–166.\cr
[17]. Hubalek Z: Coefficients of association and similarity, based on binary (presence-absence) data: An evaluation. Biol Rev 1982, 57:669–689.\cr
[18]. Jaccard P: The distribution of the flora in the alpine zone. New Phytol 1912, 11:37–50.\cr
[19]. Jackson DA, Somers KM, Harvey HH: Similarity coefficients: Measures of co-occurrence and association or simply measures of occurrence? Am Nat 1989, 133:436–453.\cr
[20]. Johnson SC: Hierarchical clustering schemes. Psychometrika 1967, 32:241–254.\cr
[21]. Lance GN, Williams WT: Computer Programs for Hierarchical Polythetic Classification (``Similarity Analyses’’). Comput J 1966, 9:60–64.\cr
[22]. Lourenco F, Lobo V, Bacao F: Binary-Based Similarity Measures for Categorical Data and Their Application in Self-Organizing Maps. 2004.\cr
[23]. Michael EL: Marine ecology and the coefficient of association: A plea in behalf of quantitative biology. J Ecol 1920, 8:54–59.\cr
[24]. Nei M, Li W-H: Mathematical model for studying genetic variation in terms of restriction endonucleases. Proc Natl Acad Sci U S A 1979, 76:5269–5273.\cr
[25]. Ojurongbe TA: Comparison of different proximity measures and classification methods for binary data. Justus Liebig University Gießen; 2012.\cr
[26]. Stiles HE: The association factor in information retrieval. J ACM 1961, 8(2):271–279.\cr
[27]. Todeschini R, Consonni V, Xiang H, Holliday J, Buscema M, Willett P: Similarity coefficients for binary chemoinformatics data: Overview and extended comparison using simulated and real data sets. J Chem Inf Model 2012, 52:2884–2901.\cr
[28]. Warrens MJ: Similarity coefficients for binary data: properties of coefficients, coefficient matrices, multi-way metrics and multivariate coefficients. Leiden University; 2008.\cr
[29]. Zhang B, Srihari SN: Binary vector dissimilarity measures for handwriting identification. In Proceedings of SPIE-IS&T Electronic Imaging Vol. 5010; 2003:28–38.\cr
[30]. Zhang B, Srihari SN: Properties of binary vector dissimilarity measures. In Proc. JCIS Int’l Conf. Computer Vision, Pattern Recognition, and Image Processing; 2003:1–4.\cr
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
%% \examples{
%% 
%% }

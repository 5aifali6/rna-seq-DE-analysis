# 🧬 RNA-Seq Differential Expression Analysis

This project performs differential gene expression (DE) analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in R, based on the [GSE60424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424) dataset from NCBI GEO.

---

## 📁 Dataset

**GSE60424**: RNA-Seq profiling of human immune cell subsets across different immune-associated diseases.

- **Organism**: *Homo sapiens*
- **Conditions compared**: LPS-stimulated vs Control
- **Source**: [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424)

---

## 🧪 Tools & Packages

- R (4.x)
- DESeq2
- pheatmap
- ggplot2
- EnhancedVolcano

---

## 📊 Results

- ✅ PCA plot (sample clustering)
- ✅ Volcano plot (log2FC vs p-value)
- ✅ Heatmap (top 30 differentially expressed genes)

> All outputs are saved in the `results/` folder.

### 🔬 Sample Plots

**Volcano Plot**  
![Volcano Plot](results/volcano_plot.png)

**PCA Plot**  
![PCA Plot](results/pca_plot.png)

**Heatmap**  
![Heatmap](results/heatmap.png)

---

## 🔁 Reproducibility

To reproduce this analysis:

```r
source("scripts/deseq2_analysis.R")

<template>
  <div class="container mx-auto p-4">
    <h1 class="text-3xl font-bold">Functions of scISD</h1>
    <p class="text-lg">
      scISD presents the analysis results in a user-friendly way for public access, which allows
      researchers to quickly gain insight into the expression of genes of interest at the
      single-cell level for various inflammatory skin diseases. It mainly includes two modules:
      Browser and Query.
    </p>

    <h1 class="text-3xl font-bold mt-4">1. Browser Module</h1>
    <p class="text-lg">
      The Browse module encompasses 34 datasets, each dataset show matched paper,disease
      type,experiment information,cell number, et al.
    </p>
    <img src="@/assets/dataset.png" />
    <p class="text-lg">
      Users can explore the cell-type-specific expression of genes of interest, identify marker
      genes for each cell type.
    </p>
    <div class="flex justify-center">
      <img src="@/assets/exp.png" />
    </div>
    <p class="text-lg">
      Furthermore, users can also view cell-type-specific gene set via select "celltype" and "Cell
      Type-specific genes". Meanwhile, significantly differentially expressed genes (DEGs) between
      disease and normal condition can also be expored via "Select Comparison Groups".
    </p>
    <img src="@/assets/deg.png" />
    <p class="text-lg">
      Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment
      analyses are provided for above significantly differentially expressed genes. A "Top 10"
      button displays a bar chart visualizing the top 10 enriched terms.
    </p>
    <img src="@/assets/go.png" />

    <h1 class="text-3xl font-bold mt-4">2. Query Module</h1>
    <p class="text-lg">
      The Query module offers two primary analysis results. Frist of all, users can input a gene of
      interest to explore its cell-type expression distribution across all datasets. A higher number
      of datasets exhibiting expression in a specific cell type suggests greater conservation of
      that gene as a marker for that cell type.
    </p>
    <img src="@/assets/cmarker.png" />
    <p class="text-lg">
      Secondly, users can input both a gene and a cell type to identify disease-conserved genes. A
      gene identified as significantly differentially expressed across multiple disease comparisons
      is considered a potential conserved gene for inflammatory skin diseases.
    </p>
    <img src="@/assets/cdeg.png" />

    <h1 class="text-3xl font-bold mt-4">3. Atlas Module</h1>
    <p class="text-lg">
      The Atlas module offers integrated data analysis results of all 34 datasets. The results are
      presented in the same manner with Browser module.
    </p>
  </div>
</template>

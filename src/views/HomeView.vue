<template>
  <!-- Main Content -->
  <main class="flex-1 bg-gray-50">
    <div class="hero bg-base-200 h-48">
      <div class="hero-content text-center">
        <div>
          <h1 class="text-7xl font-bold">scISD</h1>
          <h1 class="text-xl py-6">
            A Single-cell RNA-Seq Database for Inflammatory Skin Diseases
          </h1>
        </div>
      </div>
    </div>
    <div class="w-5/6 mx-auto py-4">
      <div>
        <h1 class="text-2xl font-bold text-center">Welcome to scISD!</h1>
        <p class="leading-relaxed text-justify p-4">
          Single-cell RNA-Seq (scRNA-seq) are an emerging method which facilitates to explore the
          comprehensive transcriptome in a single cell. To provide a comprehensive resource for
          Inflammatory Skin Diseases, we developed the scISD database, which contains 34 human
          single cell gene expression data sets collected from Gene Expression Omnibus (GEO),
          involving 1,838,748 cells from 38 cell types. We provided detailed information for gene
          expression of cells in different statuses, including boxplot of gene expression,
          differentially expression analysis and GO as well as pathway annotations for individual
          dataset. Simmilar results were also presented in a multiple-dataset integrated manner. In
          addition, we also identified common marker genes and disease-related genes via
          cross-dataset comparsion analysis.
        </p>
      </div>

      <div class="grid grid-cols-5 gap-4">
        <div class="stats shadow">
          <div class="stat">
            <div class="stat-title">Total Disease Types</div>
            <div class="stat-value">13</div>
          </div>
        </div>
        <div class="stats shadow">
          <div class="stat">
            <div class="stat-title">Total Datasets</div>
            <div class="stat-value">34</div>
          </div>
        </div>
        <div class="stats shadow">
          <div class="stat">
            <div class="stat-title">Total Samples</div>
            <div class="stat-value">483</div>
          </div>
        </div>
        <div class="stats shadow">
          <div class="stat">
            <div class="stat-title">Total Cells</div>
            <div class="stat-value">1,838,748</div>
          </div>
        </div>
        <div class="stats shadow">
          <div class="stat">
            <div class="stat-title">Total Cell Types</div>
            <div class="stat-value">38</div>
          </div>
        </div>
      </div>
      <div class="grid grid-cols-4 gap-4 py-4">
        <div class="card bg-base-100 w-96 shadow-xl">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Cell Clustering</h2>
          </div>
          <div class="card-body">
            <img src="http://scisd.databio1.com/api/images/GSE150672" class="h-64" />
          </div>
        </div>
        <div class="card bg-base-100 w-96 shadow-xl">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Gene Expression</h2>
          </div>
          <div class="card-body">
            <img src="@/assets/exp.png" class="h-64" />
          </div>
        </div>
        <div class="card bg-base-100 w-96 shadow-xl">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Function Enrichment</h2>
          </div>
          <div class="card-body">
            <img src="@/assets/goplot.png" class="h-64" />
          </div>
        </div>
        <div class="card bg-base-100 w-96 shadow-xl">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Conserved Marker</h2>
          </div>
          <div class="card-body">
            <img src="@/assets/conbar.png" class="h-64" />
          </div>
        </div>
      </div>
    </div>
  </main>
</template>
